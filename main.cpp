#include "mainwindow.h"

#include <QApplication>
#include <QtCharts>
#include <QChartView>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <math.h>
#include <QVector>
#include <QDebug>

#ifndef UseDFT

//const double pie = 3.14159;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}

#endif
