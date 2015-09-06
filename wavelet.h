#pragma once

#include <QList>
#include <QVector>

struct WaveletTransformationResult
{
    QVector<double> approx;
    QVector<double> details;
};

enum EWaveletFamily
{
    Family_Db4
};

// Forward discrete wavelet transformation
WaveletTransformationResult dwt1(const QVector<double>& signal,
                                 EWaveletFamily waveletFamily);

// Inverse discrete wavelet transformation
QVector<double> idwt1(const QVector<double>& approx, 
                      const QVector<double>& details,
                      EWaveletFamily waveletFamily,
                      int lx = -1);

QVector<double> processWavelet(const QVector<double>& data,
                               int numLevels = 7);