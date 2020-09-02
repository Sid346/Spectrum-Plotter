#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QScatterSeries>
#include <QFileDevice>
#include <QTemporaryFile>
#include <QtEndian>
#include <QVector>
#include <QSound>
#include <QMediaPlayer>
#include <QMessageBox>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    QString fileName;
    QUrl file_url;
    QVector<double> signal;
    int freq = 0;
    QMediaPlayer* player = NULL;
    ~MainWindow();
private:
    Ui::MainWindow *ui;
private slots:
    void play();
    void load();
    void on_Play_Audio_clicked();
    void on_Pause_Audio_clicked();
};
#endif // MAINWINDOW_H
