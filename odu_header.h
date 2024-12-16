#ifndef EX_32_ODU_HEADER_H
#define EX_32_ODU_HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double f(double x);
double p(double x);
double q(double x);

typedef struct Solve
{
    unsigned long long length;
    double* y;
    double* y_true;
    double lastValue;
    int n;
} Solve;

#endif //EX_32_ODU_HEADER_H
#pragma once
