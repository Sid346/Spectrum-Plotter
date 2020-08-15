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

/*void FFT(QVector<double> signal_time,
                 double* spectrum,
                 const int N
                 )
{
    int i,x,y,j,k,j1;
    float o,o1;
    complex<double> p,q,r,s,w,m,n,alpha,beta,gamma1[10],gamma2[10];
    complex<double> y1[10],b[10],c[10],d[10],gamma[10],u[10],b2[10];

    for(i = 0; i <= N/2 - 1; i++)
    {
       b[i] = signal_time[2 * i];
    }
    for(i = N/2, j = N/2 - 1; i <= N - 1,j >= 0; i++, j--)
    {
       b[i] = signal_time[i-j];
    }
    for(i = 0; i <= 1; i++)
    {
       c[i] = b[2 * i];
    }
    for(i = 2, j = 1; i <= N/2 - 1, j >= 0; i++, j--)
    {
       c[i] = b[i - j];
    }
    for(i = N / 2, j = 0; i <= N/2 + 1, j <= 1; i++, j++)
    {
        c[i] = b[i + j];
    }
    for(i=6,j=1;i<=7,j>=0;i++,j--)
    {
        c[i] = b[i-j];
    }
    for(i=0;i<=6;i=i+2)
    {
        m = c[i];
        n = c[i+1];
        alpha = m+n;
        beta = m-n;
        b[i] = alpha;
        b[i+1] = beta;
    }
    for(i=0;i<=4;i=i+4)
    {
        p = b[i];
        q = b[i+1];
        r = b[i+2];
        s = b[i+3];
        for(k=0+i,j1=0;k<=2+i,j1<=2;k=k+2,j1=j1+2)
        {
            o = cos((2*PI*j1)/4);
            o1 = sin((2*PI*j1)/4);
            w = complex<double>(o,-o1);
            gamma[k] = p+real(w)*r;
        }
        for(k=1+i,j1=1;k<=3+i,j1<=3;k=k+2,j1=j1+2)
        {
            o = cos((2*PI*j1)/4);
            o1 = sin((2*PI*j1)/4);
            w = complex<double>(o,-o1);
            gamma[k] = q+(w)*s;
        }
    }
    for(i=0;i<=7;i=i+1)
    {
          y1[i] = gamma[i];
    }

    for(k=j1=0;k<=3,j1<=3;k++,j1++)
    {
        o = cos((2*PI*j1)/8);
        o1 = sin((2*PI*j1)/8);
        w = complex<double>(o,-o1);
        gamma1[k] = y1[k]+w*y1[k+4];
    }
    for(k=j1=4;k<=7,j1<=7;k++,j1++)
    {
        o = cos((2*PI*j1)/8);
        o1 = sin((2*PI*j1)/8);
        w = complex<double>(o,-o1);
        gamma1[k] = w*y1[k]+y1[k-4];
    }
    for(i=0;i<=7;i++)
    {
        b2[i] = gamma1[i];
    }
    for(i=0;i<=7;i=i+1)
    {
        spectrum[i] = (real(b2[i]) * real(b2[i])) + (imag(b2[i]) * imag(b2[i]));
        //spectrum[i] = spectrum[i]
    }

}*/

void Swap(double *x, double *y)
{
    double temp = *x;
    *x = *y;
    *y = temp;
}

void FFT(QVector<double> signal_time, double *spectrum, long samples)
{	// Computes in-place complex-to-complex FFT
    // Reverse bits

    double* vReal = new double[samples]{ 0.0 };
    double* vImag = new double[samples]{ 0.0 };

    for(int i = 0; i < samples; i++){
        vReal[i] = signal_time[i];
    }
    int power = exp(samples);
    long j = 0;
    for (long i = 0; i < (samples - 1); i++) {
        if (i < j) {
            Swap(&vReal[i], &vReal[j]);
        }
        long k = (samples >> 1);
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    // Compute the FFT

    double c1 = -1.0;
    double c2 = 0.0;
    long l2 = 1;
    for (long l = 0; (l < power); l++) {
        long l1 = l2;
        l2 <<= 1;
        double u1 = 1.0;
        double u2 = 0.0;
        for (j = 0; j < l1; j++) {
             for (long i = j; i < samples; i += l2) {
                    long i1 = i + l1;
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
            c2 = -c2;
    }
    // Scaling for reverse transform
    for (long i = 0; i < samples; i++) {
                 vReal[i] /= samples;
                 vImag[i] /= samples;
    }
    for (long i = 0; i < samples; i++) {
            spectrum[i] = sqrt((vReal[i] * vReal[i]) + (vImag[i] * vImag[i]));
        }
}
