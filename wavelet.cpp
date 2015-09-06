#include "wavelet.h"

#include <QtAlgorithms>

namespace detail
{
    QVector<double> wavelet(EWaveletFamily family,
                            bool lowpass,
                            bool decompose)
    {
        static const double Lo_D[] = {
            -0.010597401784997,
             0.032883011666983,
             0.030841381835987,
            -0.187034811718881,
            -0.027983769416984,
             0.630880767929590,
             0.714846570552542,
             0.230377813308855
        };

        static const double Hi_D[] = {
            -0.230377813308855,
             0.714846570552542,
            -0.630880767929590,
            -0.027983769416984,
             0.187034811718881,
             0.030841381835987,
            -0.032883011666983,
            -0.010597401784997
        };

        static const double Lo_R[] = {
             0.230377813308855,
             0.714846570552542,
             0.630880767929590,
            -0.027983769416984,
            -0.187034811718881,
             0.030841381835987,
             0.032883011666983,
            -0.010597401784997
        };

        static const double Hi_R[] = {
            -0.010597401784997,
            -0.032883011666983,
             0.030841381835987,
             0.187034811718881,
            -0.027983769416984,
            -0.630880767929590,
             0.714846570552542,
            -0.230377813308855
        };

        QVector<double> ret;

        switch(family)
        {
        case Family_Db4:
            if(decompose)
            {
                if(lowpass)
                {
                    ret.resize(sizeof(Lo_D)/sizeof(double));
                    qCopy(Lo_D, Lo_D + ret.size(), ret.begin());
                }
                else
                {
                    ret.resize(sizeof(Hi_D)/sizeof(double));
                    qCopy(Hi_D, Hi_D + ret.size(), ret.begin());
                }
            }
            else
            {
                if(lowpass)
                {
                    ret.resize(sizeof(Lo_R)/sizeof(double));
                    qCopy(Lo_R, Lo_R + ret.size(), ret.begin());
                }
                else
                {
                    ret.resize(sizeof(Hi_R)/sizeof(double));
                    qCopy(Hi_R, Hi_R + ret.size(), ret.begin());
                }
            }
            break;
        }

        return ret;
    }

    QVector<double> halfPointSymmetrization(const QVector<double>& signal,
                                            int lf)
    {
        QVector<double> ret;
        ret.reserve(signal.size() + 2 * lf);
        for(int index = lf - 2; index >= 0; --index)
            ret.append(signal[index]);
        for(int index = 0; index < signal.size(); ++index)
            ret.append(signal[index]);
        for(int index = signal.size() - 1; index >= signal.size() - lf - 1; --index)
            ret.append(signal[index]);
        return ret;
    }

    QVector<double> rerange(const QVector<double>& x,
                            int first,
                            int stride,
                            int last)
    {
        Q_ASSERT(last > first);
        Q_ASSERT(x.size() >= last);

        int size = (last - first) / stride + 1;
        QVector<double> ret;
        ret.reserve(size);
        for(int index = 0; index < size; ++index)
            ret.append(x[first + index * stride]);
        return ret;
    }

    QVector<double> conv1(const QVector<double>& f,
                          const QVector<double>& g,
                          bool full)
    {
        int lf = f.size();
        int lg = g.size();
    
        if(!full)
        {
            int ly = qMax(lf - qMax(0, lg - 1), 0);
            QVector<double> y(ly, 0.0);
            for(size_t k = 0; k < ly; ++k)
            {
                for(size_t j = 0; j < lg; ++j)
                    y[k] += f[k - j + lg - 1] * g[j];
            }
            return y;
        }
        else
        {
            int ly = lf + lg - 1;
            QVector<double> y(ly, 0.0);
            for(int k = 0; k < ly; ++k)
                for(int j = 0; j < lg; ++j)
                    if(k-j >= 0 && k-j < lf)
                        y[k] += f[k - j] * g[j];
            return y;
        }
    }

    QVector<double> upsconv(const QVector<double>& x,
                            const QVector<double>& g,
                            int lx)
    {
        // Dyadic upsampling
        int p = 0;
        int rem2 = p - ((p / 2) * 2);
        int addLEN = 2 * rem2 - 1;
        int lux = 2 * x.size() + addLEN;
        QVector<double> ux(lux, 0.0);
        QVector<double>::const_iterator iter = x.begin();
        for(int index = 0 + rem2; index < ux.size(); index += 2)
            ux[index] = *iter++;

        // Convolution with recompose filter
        QVector<double> y = conv1(ux, g, true);

        // Crop boundries
        int ly = y.size();
        double d = double(ly - lx) / 2.0;
        int first = int(floor(d));
        int last = ly - int(ceil(d));

        QVector<double> z(last - first);
        qCopy(y.begin() + first, y.begin() + last, z.begin());

        return z;
    }
}

WaveletTransformationResult dwt1(const QVector<double>& signal,
                                 EWaveletFamily waveletFamily)
{
    QVector<double> Lo_D = detail::wavelet(waveletFamily, true, true);
    QVector<double> Hi_D = detail::wavelet(waveletFamily, false, true);
    Q_ASSERT(Lo_D.size() == Hi_D.size());
    int lf = Lo_D.size();

    QVector<double> y = detail::halfPointSymmetrization(signal, lf);

    int first = 1;
    int last = signal.size() + lf - 2;
    WaveletTransformationResult result;

    // Approximation
    QVector<double> z = detail::conv1(y, Lo_D, false);
    result.approx = detail::rerange(z, first, 2, last);

    // Details
    z = detail::conv1(y, Hi_D, false);
    result.details = detail::rerange(z, first, 2, last);

    return result;
}

QVector<double> idwt1(const QVector<double>& approx, 
                      const QVector<double>& details,
                      EWaveletFamily waveletFamily,
                      int lx)
{
    QVector<double> Lo_R = detail::wavelet(waveletFamily, true, false);
    QVector<double> Hi_R = detail::wavelet(waveletFamily, false, false);
    Q_ASSERT(Lo_R.size() == Hi_R.size());

    QVector<double> a = detail::upsconv(approx, Lo_R, lx);
    QVector<double> d = detail::upsconv(details, Hi_R, lx);
    Q_ASSERT(a.size() == d.size());

    QVector<double> sum(a.size());
    for(int i = 0; i < sum.size(); ++i)
        sum[i] = a[i] + d[i];   

    return sum;
}

QVector<double> processWavelet(const QVector<double>& data,
                               int numLevels)
{
    // Forward discrete wavelet transformation
    int originalSize = data.size();
    QVector<QVector<double>> details;
    QVector<double> ecgData = data;

    for(int level = 0; level < numLevels; ++level)
    {
        WaveletTransformationResult result = dwt1(ecgData, Family_Db4);
        details.append(result.details);
        ecgData = result.approx;
    }

    Q_ASSERT(details.size() == numLevels);

    // Under ecgData lies last approximation, filter them out
    for(int i = 0; i < ecgData.size(); ++i)
        ecgData[i] = 0.0;

    // Inverse discrete wavelet transformation
    QVector<double> recomp = ecgData;
    
    for(int i = 0; i < details.size(); ++i)
    {
        int index = details.size() - 1 - i;
        int recompSize = index - 1 >= 0
            ? details[index - 1].size()
            : originalSize;

        recomp = idwt1(recomp, details[index], Family_Db4, recompSize);
    }

    return recomp;
}