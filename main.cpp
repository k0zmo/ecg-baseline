#include <QVector>
#include <QList>
#include <QtAlgorithms>
#include <QFile>
#include <QTextStream>

#include "butter.h"
#include "wavelet.h"
#include "sgolay.h"

template<typename T>
void save(const T& container, const QString& filename)
{
    QFile file(filename);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream in(&file);
    for(int i = 0; i < container.size(); ++i)
        in << qSetRealNumberPrecision(20) << container[i] << "\n";
}

extern int ecg_data_int[];
extern int ecg_data_size;

int main(int argc, char *argv[])
{
    QVector<double> ecgData(ecg_data_size);
    qCopy(ecg_data_int, ecg_data_int + ecg_data_size, ecgData.begin());

    QVector<double> wv = processWavelet(ecgData, 7);
    QVector<double> bt = processButter(ecgData);
    QVector<double> sg = processSGolay(ecgData, 3, 500);

    save(ecgData, "signal_base.txt");
    save(wv, "wavelet.txt");
    save(bt, "butterworth.txt");
    save(sg, "savitzkygolay.txt");
}
