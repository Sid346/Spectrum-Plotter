#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QScatterSeries>
#include <QFileDevice>
#include <QTemporaryFile>
#include <QtEndian>
#include <QVector>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    QString fileName;
    QVector<double> signal;
    int freq = 0;
    ~MainWindow();
private:
    Ui::MainWindow *ui;
private slots:
    void play();
    void load();
};
#endif // MAINWINDOW_H
