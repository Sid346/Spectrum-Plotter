#include "test.h"
#include <QtWidgets/QApplication>
#include <QDebug>
#include<math.h>
#include<QtCharts/qchart.h>
#include<QtCharts/qchartview.h>
#include<QtCharts/qlineseries.h>
#include<QtCharts/qvalueaxis.h>


QT_CHARTS_USE_NAMESPACE

#ifndef UseDFT

const double pie = 3.14159;

int main(int argc, char* argv[])
{
	QApplication a(argc, argv);

	double freq = 400;
	

	

	test w;
	
	w.show();
	return a.exec();
}
#endif
