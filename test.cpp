#include "test.h"



extern void DFT(QVector<double> signal_time,
	double* spectrum,
	const int SIGLEN,
	double *max,
	double *min
);


test::test(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);

	int a = 10;
	connect(ui.Load, SIGNAL(clicked()), this, SLOT(play()));
	
}
extern QVector<double>* data_read(QString fileName,int* freq);

void test::play()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open Text file"), "", tr("Text Files (*.wav)"));
	int freq;
	QVector<double> signal = *data_read(fileName, &freq);
	qDebug() << "The Length of input signal is : " << signal.size() << endl;
	const int SIGLEN = signal.size();
	QChart* chart = new QChart();
	QValueAxis* axisX = new QValueAxis;
	QValueAxis* axisY = new QValueAxis;

	axisX->setRange(0, freq/2);
	
	chart->addAxis(axisX, Qt::AlignBottom);
	chart->addAxis(axisY, Qt::AlignLeft);
	QLineSeries* series0 = new QLineSeries();
	series0->setName("scatter1");
	int N = 1024;
	double* spectrum = new double[N / 2]{ 0.0 };
	int i = 0;
	double maximum = -DBL_MAX,minimum = DBL_MAX,max, min;
	double* sum = new double[N / 2]{ 0.0 };
	int  index = 0;

	while ( i < SIGLEN-N) {
		QVector<double> window(signal.begin()+i, signal.begin() + i+ N);
		max = -DBL_MAX, min = DBL_MAX;

		DFT(window, spectrum, N, &max,&min);
		for (int n = 0; n < N / 2; n++) {
			sum[n] +=spectrum[n]/SIGLEN;
		}
		i += N;
	}
	for (int n = 0; n < N / 2; n++) {
		if (sum[n] > maximum) {
			maximum = sum[n];
			index = n;
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
	
	chart->addSeries(series0);
	
	series0->attachAxis(axisX);
	series0->attachAxis(axisY);
	chart->setTitle("Spectrum Plotter");
	chart->setDropShadowEnabled(false);
	QChartView* view = new QChartView(chart);
	view->setRenderHint(QPainter::Antialiasing);
	ui.chartview->setChart(chart);

}
