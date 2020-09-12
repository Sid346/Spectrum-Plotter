#include "mainwindow.h"
#include "ui_mainwindow.h"

extern void DFT(QVector<double> signal_time,
    double* spectrum,
    const int SIGLEN
);

/*extern void FFT(QVector<double> signal_time,
    double* spectrum,
    const int SIGLEN
);*/
extern void FFT(QVector<double> signal_time, double *spectrum, long samples);

extern void hamming(QVector<double>& signal_time, const int N);
extern QVector<double>* LMS(QVector<double> x, QVector<double> d, float delta, int L);

//QVector<double> signal;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->progressBar->setValue(0);
    connect(ui->Load, SIGNAL(clicked()), this, SLOT(load()));
    connect(ui->Calculate, SIGNAL(clicked()), this, SLOT(play()));
}

extern QVector<double>* data_read(QString fileName,int* freq);
extern void data_write(QVector<double> data, QString filename);

void MainWindow::load() {

    MainWindow::fileName = QFileDialog::getOpenFileName(this,
        tr("Open Text file"), "", tr("Text Files (*.wav)"));

    player = new QMediaPlayer;
    MainWindow::file_url = QUrl::fromLocalFile(MainWindow::fileName);
    player->setMedia(file_url);

    ui->progressBar->setValue(0);
    MainWindow::signal = *data_read(MainWindow::fileName, &(MainWindow::freq));
    ui->progressBar->setValue(20);
    const int SIGLEN = MainWindow::signal.size();
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
        ui->progressBar->setValue(i*80/SIGLEN);
        series0->append(i, signal[i]);
        i++;
    }
    ui->progressBar->setValue(80);
    axisY->setRange(-1, 1);
    chart->addSeries(series0);
    chart->legend()->setVisible(false);
    series0->attachAxis(axisX);
    ui->progressBar->setValue(90);
    series0->attachAxis(axisY);

    ui->Signalview->setChart(chart);
    ui->progressBar->setValue(100);
}

void MainWindow::play()
{
    if(ui->radioButton->isChecked()){
        ui->progressBar->setValue(0);
        const int SIGLEN = MainWindow::signal.size();
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
            ui->progressBar->setValue(i*100/(SIGLEN - N));
            QVector<double> window;
            window.resize(N);
            for(int j = 0; j < N; j++){
                window[j] = *(MainWindow::signal.begin() + j + i);
            }
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
        ui->progressBar->setValue(100);
        ui->Spectrumview->setChart(chart);

    }

    if(ui->radioButton_2->isChecked()){
        ui->progressBar->setValue(0);
        const int SIGLEN = MainWindow::signal.size();
        QChart* chart = new QChart();
        QValueAxis* axisX = new QValueAxis;
        QValueAxis* axisY = new QValueAxis;

        axisX->setRange(0, freq / 2);

        chart->addAxis(axisX, Qt::AlignBottom);
        chart->addAxis(axisY, Qt::AlignLeft);
        QLineSeries* series0 = new QLineSeries();
        int N = 1024;
        double* spectrum = new double[N]{ 0.0 };
        int i = 0;
        double maximum = -DBL_MAX, minimum = DBL_MAX;
        double* sum = new double[N / 2]{ 0.0 };

        while ( i < SIGLEN-N) {
            ui->progressBar->setValue(i*100/(SIGLEN - N));
            QVector<double> window;
            window.resize(N);
            for(int j = 0; j < N; j++){
                window[j] = *(MainWindow::signal.begin() + j + i);
            }
            hamming(window, N);
            FFT(window, spectrum, N);

            for (int n = 0; n < N / 2; n++) {
                sum[n] += spectrum[n] / SIGLEN;
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
        ui->progressBar->setValue(100);
        ui->Spectrumview->setChart(chart);
    }

}

void MainWindow::on_Play_Audio_clicked()
{
    //QSound audio(MainWindow::fileName);
    //audio.play();
    if(player == NULL){
        QMessageBox::about(this, "Error", "Please select audio file first");
    }
    else
        player->play();
}

void MainWindow::on_Pause_Audio_clicked()
{
    player->pause();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_OK_clicked()
{
    if(ui->comboBox->currentText() == "LMS"){
        //int fs = MainWindow::freq;
        int L = 16;
        int M = MainWindow::signal.size();
        int delay = L;
        float delta = 0.005;
        QVector<double> zd;
        zd.resize(M);
        for(int i = 0; i < delay - 1; i++){
            zd[i] = 0;
        }
        for(int j = delay - 1; j < M; j++){
            zd[j] = *(MainWindow::signal.begin() + j - delay + 1);
        }

        ui->progressBar->setValue(0);

        QVector<double> y;
        y.resize(M);
        y = *LMS(zd, MainWindow::signal, delta, L);

        data_write(y, "LMS_output.wav");

        player->setMedia(QUrl("LMS_output.wav"));

        ui->progressBar->setValue(20);

        const int SIGLEN = y.size();
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
            ui->progressBar->setValue(i*80/SIGLEN);
            series0->append(i, y[i]);
            i++;
        }
        ui->progressBar->setValue(80);
        axisY->setRange(-1, 1);
        chart->addSeries(series0);
        chart->legend()->setVisible(false);
        series0->attachAxis(axisX);
        ui->progressBar->setValue(90);
        series0->attachAxis(axisY);

        ui->Signalview->setChart(chart);
        ui->progressBar->setValue(100);
    }
}
