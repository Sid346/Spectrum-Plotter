#include "test.h"



extern void DFT(QVector<double> signal_time,
	double* spectrum,
	const int SIGLEN
);

extern void hamming(QVector<double>& signal_time, const int N);

test::test(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	connect(ui.Load, SIGNAL(clicked()), this, SLOT(load()));
	connect(ui.Calculate, SIGNAL(clicked()), this, SLOT(play()));
	
}
extern QVector<double>* data_read(QString fileName,int* freq);

void test::load() {
	test::fileName = QFileDialog::getOpenFileName(this,
		tr("Open Text file"), "", tr("Text Files (*.wav)"));
	
	test::signal = *data_read(test::fileName, &(test::freq));
	const int SIGLEN = test::signal.size();
	QChart* chart = new QChart();
	QValueAxis* axisX = new QValueAxis;
	QValueAxis* axisY = new QValueAxis;
	axisX->setRange(0, SIGLEN);
	axisX->setMinorTickCount(1);
	axisX->setTickCount(10);
	chart->addAxis(axisX, Qt::AlignBottom);
	chart->addAxis(axisY, Qt::AlignLeft);
	QLineSeries* series0 = new QLineSeries();
	int i = 0;
	while (i < SIGLEN) {
		series0->append(i, signal[i]);
		i++;
	}
	axisY->setRange(-1, 1);
	chart->addSeries(series0);
	chart->legend()->setVisible(false);
	series0->attachAxis(axisX);
	
	series0->attachAxis(axisY);

	ui.Signalview->setChart(chart);

}

void test::play()
{
	const int SIGLEN = test::signal.size();
	QChart* chart = new QChart();
	QValueAxis* axisX = new QValueAxis;
	QValueAxis* axisY = new QValueAxis;

	axisX->setRange(0, freq/2);
	
	chart->addAxis(axisX, Qt::AlignBottom);
	chart->addAxis(axisY, Qt::AlignLeft);
	QLineSeries* series0 = new QLineSeries();
	int N = 1024;
	double* spectrum = new double[N / 2]{ 0.0 };
	int i = 0;
	double maximum = -DBL_MAX, minimum = DBL_MAX;
	double* sum = new double[N / 2]{ 0.0 };

	while ( i < SIGLEN-N) {
		QVector<double> window(signal.begin()+i, signal.begin() + i+ N);
		hamming(window, N);
		DFT(window, spectrum, N);
		for (int n = 0; n < N / 2; n++) {
			sum[n] += spectrum[n]/SIGLEN;
		}
		i += N;
	}
	for (int n = 0; n < N / 2; n++) {
		if (sum[n] > maximum) {
			maximum = sum[n];
		}

		if (sum[n] < minimum) {
			minimum = sum[n];
		}
	}
	for (int n = 0; n < N / 2; n++)
	{
		series0->append(((n) * (freq / 2) / (N / 2 - 1)), 10*log10(sum[n]));
	}
	axisY->setRange(10*log10(minimum), 10*log10(maximum));
	axisX->setMinorTickCount(1);
	axisX->setTickCount(10);
	chart->addSeries(series0);
	chart->legend()->setVisible(false);
	series0->attachAxis(axisX);
	series0->attachAxis(axisY);

	ui.Spectrumview->setChart(chart);

}
