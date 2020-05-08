#include "test.h"

const int SIGLEN = 10000;

extern void DFT(const double* signal_time,
	double* spectrum,
	const int SIGLEN
);


test::test(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	int a = 10;
	connect(ui.Load, SIGNAL(clicked()), this, SLOT(play()));
	
}

void test::play()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open Text file"), "", tr("Text Files (*.txt)"));
	QFile file(fileName);
	file.open(QIODevice::ReadOnly);
	QTextStream in(&file);
	double spectrum[SIGLEN / 2] = { 0.0 };
	double signal[SIGLEN] = { 0.0 };
	int i = 0;
	while (!in.atEnd())
	{
			
			QString line = in.readLine();
			QTextStream lineStream(&line);
			lineStream >> signal[i];
			qDebug() << signal[i] << endl;
			i++;

	}
	file.close();

	QLineSeries* series0 = new QLineSeries();
	series0->setName("scatter1");

	DFT(signal, spectrum, SIGLEN);
	for (int i = 0; i < SIGLEN / 2; i++) {
		series0->append(i * 22050 / SIGLEN * 2, spectrum[i] / 1000);
	}
	QChart* chart = new QChart();
	QValueAxis* axisX = new QValueAxis;
	QValueAxis* axisY = new QValueAxis;
	axisX->setRange(0, 5000);
	axisY->setRange(0, 1);
	chart->addSeries(series0);
	chart->addAxis(axisX, Qt::AlignBottom);
	series0->attachAxis(axisX);
	chart->addAxis(axisY, Qt::AlignLeft);
	series0->attachAxis(axisY);
	chart->setTitle("Spectrum Plotter");
	chart->setDropShadowEnabled(false);
	QChartView* view = new QChartView(chart);
	view->setRenderHint(QPainter::Antialiasing);
	ui.chartview->setChart(chart);

}
