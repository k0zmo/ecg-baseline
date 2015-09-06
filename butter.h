#pragma once

#include <QVector>
#include <complex>

typedef std::complex<double> complex;

// Second-order (biquadratic) IIR filtering
QVector<double> sosfilt(const QVector<double>& sosMatrix,
                        double gain,
                        const QVector<double>& x);
// Second-order (biquadratic) IIR Zero-phase digital filtering
QVector<double> sosfiltfilt(const QVector<double>& sosMatrix,
                            float gain,
                            const QVector<complex>& poles,
                            const QVector<double>& x);
QVector<double> processButter(const QVector<double>& signal);
