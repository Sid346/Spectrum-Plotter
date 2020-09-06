#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_test.h"
#include <QtCharts/qscatterseries.h>


class test : public QMainWindow
{
	Q_OBJECT

public:
	test(QWidget *parent = Q_NULLPTR);
	QString fileName;
	QVector<double> signal;
	int freq = 0;
private:
	Ui::testClass ui;
	
private slots:
	void play();
	void load();
	
};
