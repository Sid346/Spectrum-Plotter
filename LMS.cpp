#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include<iostream>
#include<conio.h>
#include<dos.h>
#include<dos.h>
#include<complex>
#include <QPoint>
#include <QtCore>

#pragma warning(disable:4996)


QVector<double>* LMS(QVector<double> x, QVector<double> d, float delta, int L){
    int N = x.size();
    QVector<double> b, y, x1, e;
    QVector<double>* y1 = new QVector<double>;
    b.resize(L);
    x1.resize(L);
    e.resize(N);
    for(int i = 0; i < L; i++){
        b[i] = 0;
        x1[i] = 0;
    }
    y.resize(N);
    for(int i = 0; i < N; i++){
        y[i] = 0;
        e[i] = 0;
    }

    for(int n = L - 1; n < N; n++){
        int j = 0;
        for(int i = n; i >= n - L + 1; i--){
            x1[j] = x[i];
            j++;
        }

        for(int i = 0; i < L; i++){
            y[n] += b[i] * x1[i];
        }

        e[n] = d[n] - y[n];

        for(int i = 0; i < L; i++){
            b[i] += delta * e[n] * x1[i];
        }
    }
    for(int i = 0; i < N; i++){
        *y1 << y[i];
    }
    return y1;
}

