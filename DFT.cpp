#include "test.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include<iostream>




#pragma warning(disable:4996)

#define PIE         3.14159

using namespace std;


void DFT(QVector<double> signal_time,
                 double* spectrum,
                 const int SIGLEN
                 )
{
 
    double* out_real = new double[SIGLEN / 2];
    double* out_imag = new double[SIGLEN / 2];
    for (int i = 0; i < SIGLEN / 2; i++) {
        out_real[i] = 0.0;
        out_imag[i] = 0.0;
    }
    for (int i = 0; i < SIGLEN /2; i++) {
        for (int j = 0; j < SIGLEN; j++) {

                out_real[i] += signal_time[j] * cos(2*PIE*i*j/SIGLEN);
                out_imag[i] += signal_time[j] * -sin(2 * PIE * i * j /SIGLEN);
            }
    }
    for (int i = 0; i < SIGLEN / 2; i++) {
        spectrum[i] = sqrt(pow(out_real[i], 2) + pow(out_imag[i], 2));
    }

}

