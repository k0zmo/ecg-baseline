#include "sgolay.h"

#include <Eigen/LU>

Eigen::MatrixXd invert(const Eigen::MatrixXd& A)
{
    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(A.rows(), A.cols());
    Eigen::PartialPivLU<Eigen::MatrixXd> LU = A.lu();
    return LU.solve(B);
}

Eigen::VectorXd sgolayfiltCoeff(const Eigen::VectorXd& x,
                                const int deg)
{
    int rows(x.size());
    int cols(deg + 1);
    Eigen::MatrixXd A(rows, cols);

    // generate input matrix for least squares fit
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            A(i,j) = pow(double(i), double(j));

    Eigen::MatrixXd AT(A.transpose());
    Eigen::MatrixXd c(invert(AT*A) * (AT * x));

    Eigen::VectorXd res(rows);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < rows; ++i)
    {
        res(i) = c(0, 0);
        for(int j = 1; j < cols; ++j)
            res(i) += c(j, 0) * A(i,j);
    }

    return res;
}

Eigen::VectorXd sgolayfilt(const Eigen::VectorXd& x,
                           const int deg,
                           const int frame)
{
    const int window = 2 * frame + 1;
    const int endidx = x.size() - 1;

    Eigen::VectorXd res(Eigen::VectorXd::Zero(x.size()));
    if(frame < 1 || deg < 1 || x.size() < (2*frame + 2))
        return res;

    // Borders
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < frame; ++i)
    {
        Eigen::VectorXd b1(Eigen::VectorXd::Zero(window));
        b1(i) = 1.0;

        const Eigen::VectorXd c1(sgolayfiltCoeff(b1, deg));
        for(int j = 0; j < window; ++j)
        {
            res(i)          += c1(j) * x(j);
            res(endidx - i) += c1(j) * x(endidx - j);
        }
    }

    // Rest of signal
    Eigen::VectorXd b2(Eigen::VectorXd::Zero(window));
    b2(frame) = 1.0;
    const Eigen::VectorXd c2(sgolayfiltCoeff(b2, deg));

    #pragma omp parallel for schedule(static)
    for(int i = 0; i <= (x.size() - window); ++i)
        for(int j = 0; j < window; ++j)
            res(i + frame) += c2(j) * x(i + j);
    return res;
}

QVector<double> processSGolay(const QVector<double>& signal,
                              const int deg,
                              const int frame,
                              QVector<double>* baselineMode)
{
    Eigen::VectorXd v(signal.size());
    qCopy(signal.begin(), signal.end(), &v(0));

    Eigen::VectorXd inter = sgolayfilt(v, deg, frame);
    v = v - inter;

    if(baselineMode != 0)
    {
        baselineMode->resize(v.size());
        qCopy(&inter(0), &inter(inter.size()-1)+1, baselineMode->begin());
    }

    QVector<double> ret(v.size());
    qCopy(&v(0), &v(v.size()-1)+1, ret.begin());

    return ret;
}
