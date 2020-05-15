#include "test.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include<iostream>




#pragma warning(disable:4996)

#define PI        3.14159

using namespace std;


void DFT(QVector<double> signal_time,
                 double* spectrum,
                 const int N,
                 double *max,
                 double *min
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
        if (spectrum[i] > *max)
            *max = spectrum[i];
        if (spectrum[i] < *min)
            *min = spectrum[i];
    }

}
//void fft(QVector<double> signal_time,
//    double* spectrum,
//    const int SIGLEN)
//
//inline void fft(long m, double* x, double* y)
//{
//    long n, i, i1, j, k, i2, l, l1, l2;
//    double c1, c2, tx, ty, t1, t2, u1, u2, z;
//
//    /* Calculate the number of points */
//    n = 1;
//    for (i = 0; i < m; i++)
//        n *= 2;
//
//    /* Do the bit reversal */
//    i2 = n >> 1;
//    j = 0;
//    for (i = 0; i < n - 1; i++) {
//        if (i < j) {
//            tx = x[i];
//            ty = y[i];
//            x[i] = x[j];
//            y[i] = y[j];
//            x[j] = tx;
//            y[j] = ty;
//        }
//        k = i2;
//        while (k <= j) {
//            j -= k;
//            k >>= 1;
//        }
//        j += k;
//    }
//
//    /* Compute the FFT */
//    c1 = -1.0;
//    c2 = 0.0;
//    l2 = 1;
//    for (l = 0; l < m; l++) {
//        l1 = l2;
//        l2 <<= 1;
//        u1 = 1.0;
//        u2 = 0.0;
//        for (j = 0; j < l1; j++) {
//            for (i = j; i < n; i += l2) {
//                i1 = i + l1;
//                t1 = u1 * x[i1] - u2 * y[i1];
//                t2 = u1 * y[i1] + u2 * x[i1];
//                x[i1] = x[i] - t1;
//                y[i1] = y[i] - t2;
//                x[i] += t1;
//                y[i] += t2;
//            }
//            z = u1 * c1 - u2 * c2;
//            u2 = u1 * c2 + u2 * c1;
//            u1 = z;
//        }
//        c2 = sqrt((1.0 - c1) / 2.0);
//        
//        c2 = -c2;
//        c1 = sqrt((1.0 + c1) / 2.0);
//    }
//
//    /* Scaling for forward transform */
// 
//        for (i = 0; i < n; i++) {
//            x[i] /= n;
//            y[i] /= n;
//        }
//  
//}
