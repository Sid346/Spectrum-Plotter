#include "mainwindow.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include<iostream>
#include<conio.h>
#include<dos.h>
#include<complex>
#include <QPoint>
#include <QtCore>

#pragma warning(disable:4996)

#define PI        3.14159

using namespace std;

void hamming(QVector<double> &signal_time, const int N) {
    for (int i = 0; i < N; i++)
        signal_time[i] *= (0.54 - 0.46 * cos(2*PI * i/N));
}


void DFT(QVector<double> signal_time,
                 double* spectrum,
                 const int N
                 )
{

    double* out_real = new double[N / 2]{ 0 };
    double* out_imag = new double[N / 2]{ 0 };

    for (int i = 0; i < N /2; i++) {
        for (int j = 0; j < N; j++) {
                out_real[i] += signal_time[j] * cos(2*PI*i*j/N);
                out_imag[i] += signal_time[j] * -sin(2 * PI * i * j /N);
            }
        spectrum[i] = out_real[i]* out_real[i] + out_imag[i]* out_imag[i];
    }

}


void Swap(double *x, double *y)
{
    double temp = *x;
    *x = *y;
    *y = temp;
}
uint8_t exponent(uint16_t value){
    uint8_t result = 0;
        while (((value >> result) & 1) != 1) result++;
        return(result);
}


void FFT(double *vReal, double *vImag ,long samples, uint8_t dir)
{	// Computes in-place complex-to-complex FFT
    // Reverse bits


    uint8_t power = exponent(samples);

    uint16_t j = 0;
    for (uint16_t i = 0; i < (samples - 1); i++) {
        if (i < j) {
            Swap(&vReal[i], &vReal[j]);
            if(dir == FFT_REVERSE)
                            Swap(&vImag[i], &vImag[j]);
        }
        uint16_t k = (samples >> 1);
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    // Compute the FFT
    double c1 = -1.0;
    double c2 = 0.0;
    uint16_t l2 = 1;
    for (uint8_t l = 0; (l < power); l++) {
        uint16_t l1 = l2;
        l2 <<= 1;
        double u1 = 1.0;
        double u2 = 0.0;
        for (j = 0; j < l1; j++) {
           for (uint16_t i = j; i < samples; i += l2) {
                uint16_t i1 = i + l1;
                double t1 = u1 * vReal[i1] - u2 * vImag[i1];
                double t2 = u1 * vImag[i1] + u2 * vReal[i1];
                vReal[i1] = vReal[i] - t1;
                vImag[i1] = vImag[i] - t2;
                vReal[i] += t1;
                vImag[i] += t2;
            }
            double z = ((u1 * c1) - (u2 * c2));
            u2 = ((u1 * c2) + (u2 * c1));
            u1 = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        c1 = sqrt((1.0 + c1) / 2.0);
        if (dir == FFT_FORWARD) {
            c2 = -c2;
        }

    }
    if (dir != FFT_FORWARD) {
        for (uint16_t i = 0; i < samples; i++) {
             vReal[i] /= samples;
             vImag[i] /= samples;
        }
    }

}
