#pragma once

#include <QVector>

QVector<double> processSGolay(const QVector<double>& ecgData,
                              const int deg = 3,
                              const int frame = 500,
                              QVector<double>* baselineModel = 0);
