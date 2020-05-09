#include "test.h"



extern void DFT(QVector<double> signal_time,
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
extern QVector<double>* data_read(QString fileName);

void test::play()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open Text file"), "", tr("Text Files (*.wav)"));

	QVector<double> signal = *data_read(fileName);
	qDebug() << "The Length of input signal is : " << signal.size() << endl;
	const int SIGLEN = signal.size();
	QLineSeries* series0 = new QLineSeries();
	series0->setName("scatter1");
	double *spectrum = new double[SIGLEN/2];
	DFT(signal, spectrum, SIGLEN);
	for (int i = 0; i < SIGLEN / 2; i++) {
		series0->append(i * 22050 / (SIGLEN / 2), spectrum[i] / 100);
	}
	QChart* chart = new QChart();
	QValueAxis* axisX = new QValueAxis;
	QValueAxis* axisY = new QValueAxis;
	axisX->setRange(0, 22050);
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
