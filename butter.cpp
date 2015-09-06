#include "butter.h"

#include <QtAlgorithms>
#include <cmath>
#include <limits>

namespace detail
{
double abs(const complex& c)
{
    return std::sqrt(c.real()*c.real() + c.imag()*c.imag());
}

double angle(const complex& c)
{
    return std::atan2(c.imag(), c.real());
}
}

QVector<double> sosfilt(const QVector<double>& sosMatrix,
                        double gain,
                        const QVector<double>& x)
{
    Q_ASSERT(!(sosMatrix.size() % 5));

    int length = x.size();
    QVector<double> y(length, 0);
    QVector<double> tmp = x;
    int sosLength = sosMatrix.size() / 5;

    for(int k = 0; k < sosLength; ++k)
    {
        double v1 = 0.0f, v2 = 0.0f;

        double b0 = sosMatrix[k * 5 + 0];
        double b1 = sosMatrix[k * 5 + 1];
        double b2 = sosMatrix[k * 5 + 2];
        double a1 = sosMatrix[k * 5 + 3];
        double a2 = sosMatrix[k * 5 + 4];

        for(int n = 0; n < length; ++n)
        {
            // Assuming a0 = 1:
            // y[n] + a1*y[n-1] + a2*y[n-2] = b0*x[n] + b1*x[n-1] + b2*x[n-2]
            double v0 = tmp[n] - a1*v1 - a2*v2;
            y[n] = b0*v0 + b1*v1 + b2*v2;
            v2 = v1;
            v1 = v0;
        }

        tmp = y;
    }

    if(fabsf(gain - 1.0f) > std::numeric_limits<double>::epsilon())
    {
        for(int n = 0; n < length; ++n)
            y[n] *= gain;
    }

    return y;
}

QVector<double> sosfiltfilt(const QVector<double>& sosMatrix,
                            double gain,
                            const QVector<complex>& poles,
                            const QVector<double>& x)
{
    int transient = 0;
    static const double PI = 3.141592653589793f;
    int l = x.size();

    // Approximate duration of startup transient
    for(int i = 0; i < poles.size(); ++i)
    {
        double t = PI / (detail::abs(poles[i]) * detail::angle(poles[i]));
        transient = int(ceilf(qMax(double(transient), t)));
    }
    transient = qMin(transient, l - 1);

    // Pad with reflected data to reduce startup transients
    double pre = x[0];
    double suf = x[l-1];

    QVector<double> xx;
    xx.reserve(transient*2 + l);
    for(int index = transient; index >= 1; --index)
        xx.append(2 * pre - x[index]);
    for(int index = 0; index < l; ++index)
        xx.append(x[index]);
    for(int index = l-2; index >= l-transient-1; --index)
        xx.append(2 * suf - x[index]);

    // Forward filter data
    QVector<double> y = sosfilt(sosMatrix, 1.0f, xx);

    // Reverse result
    int ll = y.size();
    QVector<double> yflip(ll);
    for(int i = 0; i < ll; ++i)
        yflip[i] = y[ll - i - 1];

    // Reverse filter data
    yflip = sosfilt(sosMatrix, 1.0f, yflip);

    // Reverse result
    for(int i = 0; i < ll; ++i)
        y[i] = yflip[ll - i - 1];

    // Remove transients
    QVector<double> tmp = y;
    y.resize(l);
    qCopy(tmp.begin() + transient, tmp.begin() + tmp.size() - transient, y.begin());


    if(fabsf(gain - 1.0f) > std::numeric_limits<double>::epsilon())
    {
        for(int n = 0; n < l; ++n)
            y[n] *= gain;
    }

    return y;
}

double sos_butter_matrix[] = {
    0.993848328562109, -1.987696657124219, 0.993848328562109, -1.987658813704708, 0.987734500543730,
    0.991435680867689, -1.982871361735378, 0.991435680867689, -1.982833610182527, 0.982909113285588,
    0.778659906650650, -1.557319813301300, 0.778659906650650, -1.997197410205620, 0.997273460258999,
    1.000000000000000, -2.000000000000000, 1.000000000000000, -1.984493614355032, 0.984569180668372,
    1.000000000000000, -2.000000000000000, 1.000000000000000, -1.992031885495420, 0.992107738853956
};
const int sos_butter_matrix_size = sizeof(sos_butter_matrix) / sizeof(double);
double sos_butter_gain = 1.267517578260692;

complex sos_butter_poles[] = {
    complex(0.993829406852354, +0.006132749728139),
    complex(0.993829406852354, -0.006132749728139),
    complex(0.991416805091264, +0.001353465263365),
    complex(0.991416805091264, -0.001353465263365),
    complex(0.998598705102810, +0.008607347209794),
    complex(0.998598705102810, -0.008607347209794),
    complex(0.992246807177516, +0.003931197577972),
    complex(0.992246807177516, -0.003931197577972),
    complex(0.996015942747710, +0.007744717318717),
    complex(0.996015942747710, -0.007744717318717)
};
const int sos_butter_poles_size = sizeof(sos_butter_poles) / sizeof(complex);

QVector<double> processButter(const QVector<double>& signal)
{
    double sosGain = sos_butter_gain;
    QVector<double> sosMatrix(sos_butter_matrix_size);
    qCopy(sos_butter_matrix, sos_butter_matrix + sos_butter_matrix_size, sosMatrix.begin());

    QVector<complex> poles(sos_butter_poles_size);
    qCopy(sos_butter_poles, sos_butter_poles + sos_butter_poles_size, poles.begin());

    return sosfiltfilt(sosMatrix, sosGain, poles, signal);
}
