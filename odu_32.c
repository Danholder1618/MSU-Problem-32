#include "odu_header.h"

double F1(double y2)
{
    return y2;
}

double F2(double x, double y1, double y2)
{
    return -p(x) * y2 - q(x) * y1 + f(x);
}

double TrueAnswer(double x)
{
    return x / (x * x + 1);
}

void WriteTrueAnswer(double a, double h, Solve* solve)
{
    double x = a;
    solve->y_true = (double*)malloc(solve->length * sizeof(double));
    for (int i = 0; i < solve->length; i++) {
        solve->y_true[i] = TrueAnswer(x);
        x += h;
    }
}

void RungeKutta(double* y1, double* y2, double t, double h, int i) {
    double k1_y1 = h * F1(y2[i - 1]);
    double k1_y2 = h * F2(t, y1[i - 1], y2[i - 1]);

    double k2_y1 = h * F1(y2[i - 1] + k1_y2 / 4);
    double k2_y2 = h * F2(t + h / 4, y1[i - 1] + k1_y1 / 4, y2[i - 1] + k1_y2 / 4);

    double k3_y1 = h * F1(y2[i - 1] + (3 * k1_y2 + 9 * k2_y2) / 32);
    double k3_y2 = h * F2(t + 3 * h / 8, y1[i - 1] + (3 * k1_y1 + 9 * k2_y1) / 32, y2[i - 1] + (3 * k1_y2 + 9 * k2_y2) / 32);

    double k4_y1 = h * F1(y2[i - 1] + (1932 * k1_y2 - 7200 * k2_y2 + 7296 * k3_y2) / 2197);
    double k4_y2 = h * F2(t + 12 * h / 13, y1[i - 1] + (1932 * k1_y1 - 7200 * k2_y1 + 7296 * k3_y1) / 2197,
        y2[i - 1] + (1932 * k1_y2 - 7200 * k2_y2 + 7296 * k3_y2) / 2197);

    double k5_y1 = h * F1(y2[i - 1] + (439 * k1_y2 / 216 - 8 * k2_y2 + 3680 * k3_y2 / 513 - 845 * k4_y2 / 4104));
    double k5_y2 = h * F2(t + h, y1[i - 1] + (439 * k1_y1 / 216 - 8 * k2_y1 + 3680 * k3_y1 / 513 - 845 * k4_y1 / 4104),
        y2[i - 1] + (439 * k1_y2 / 216 - 8 * k2_y2 + 3680 * k3_y2 / 513 - 845 * k4_y2 / 4104));

    double k6_y1 = h * F1(y2[i - 1] - (8 * k1_y2 / 27 + 2 * k2_y2 - 3544 * k3_y2 / 2565 + 1859 * k4_y2 / 4104 - 11 * k5_y2 / 40));
    double k6_y2 = h * F2(t + h / 2, y1[i - 1] - (8 * k1_y1 / 27 + 2 * k2_y1 - 3544 * k3_y1 / 2565 + 1859 * k4_y1 / 4104 - 11 * k5_y1 / 40),
        y2[i - 1] - (8 * k1_y2 / 27 + 2 * k2_y2 - 3544 * k3_y2 / 2565 + 1859 * k4_y2 / 4104 - 11 * k5_y2 / 40));

    y1[i] = y1[i - 1] + (16 * k1_y1 / 135 + 6656 * k3_y1 / 12825 + 28561 * k4_y1 / 56430 - 9 * k5_y1 / 50 + 2 * k6_y1 / 55);
    y2[i] = y2[i - 1] + (16 * k1_y2 / 135 + 6656 * k3_y2 / 12825 + 28561 * k4_y2 / 56430 - 9 * k5_y2 / 50 + 2 * k6_y2 / 55);
}


Solve* AdamsMoulton(double y1_0, double y2_0, double a, double b, double h, Solve* solve, int needSolveWithPoints) {
    unsigned long long pointsCount = fabs(b - a) / h + 1;
    double* y1 = (double*)malloc(pointsCount * sizeof(double));
    double* y2 = (double*)malloc(pointsCount * sizeof(double));

    y1[0] = y1_0;
    y2[0] = y2_0;
    double t = a;

    for (int i = 1; i < 6; i++) {
        RungeKutta(y1, y2, t, h, i);
        t += h;
    }

    for (int i = 5; i < pointsCount - 1; i++) {
        double predictor_y1 = y1[i] + h / 1440.0 * (4277 * F1(y2[i]) - 7923 * F1(y2[i - 1]) +
            9982 * F1(y2[i - 2]) - 7298 * F1(y2[i - 3]) +
            2877 * F1(y2[i - 4]) - 475 * F1(y2[i - 5]));

        double predictor_y2 = y2[i] + h / 1440.0 * (4277 * F2(t, y1[i], y2[i]) - 7923 * F2(t - h, y1[i - 1], y2[i - 1]) +
            9982 * F2(t - 2 * h, y1[i - 2], y2[i - 2]) - 7298 * F2(t - 3 * h, y1[i - 3], y2[i - 3]) +
            2877 * F2(t - 4 * h, y1[i - 4], y2[i - 4]) - 475 * F2(t - 5 * h, y1[i - 5], y2[i - 5]));

        double corrector_y1, corrector_y2;
        for (int iter = 0; iter < 3; iter++) {
            corrector_y1 = y1[i] + h / 1440.0 * (475 * F1(predictor_y2) + 1427 * F1(y2[i]) -
                798 * F1(y2[i - 1]) + 482 * F1(y2[i - 2]) -
                173 * F1(y2[i - 3]) + 27 * F1(y2[i - 4]));

            corrector_y2 = y2[i] + h / 1440.0 * (475 * F2(t + h, predictor_y1, predictor_y2) + 1427 * F2(t, y1[i], y2[i]) -
                798 * F2(t - h, y1[i - 1], y2[i - 1]) + 482 * F2(t - 2 * h, y1[i - 2], y2[i - 2]) -
                173 * F2(t - 3 * h, y1[i - 3], y2[i - 3]) + 27 * F2(t - 4 * h, y1[i - 4], y2[i - 4]));

            predictor_y1 = corrector_y1;
            predictor_y2 = corrector_y2;
        }

        y1[i + 1] = corrector_y1;
        y2[i + 1] = corrector_y2;
        t += h;
    }

    if (needSolveWithPoints == 1) {
        solve->y = (double*)malloc(pointsCount * sizeof(double));
        for (int i = 0; i < pointsCount; i++) {
            solve->y[i] = y1[i];
        }
        solve->lastValue = y1[pointsCount - 1];
        solve->length = pointsCount;
        free(y2);
        free(y1);
    }
    else {
        solve->y = NULL;
        solve->lastValue = y1[pointsCount - 1];
        solve->length = pointsCount;
        free(y2);
        free(y1);
    }
    return solve;
}

Solve* ShootingMethod(double y_in_a, double y_in_b, double a, double b, double h, double e)
{
    Solve* y = (Solve*)malloc(sizeof(Solve));
    y->n = 0;
    Solve* temp = (Solve*)malloc(sizeof(Solve));

    double parameterLow = -10.0;
    double parameterHigh = 10.0;

    temp = AdamsMoulton(y_in_a, parameterLow, a, b, h, temp, 0);
    double lowValue = temp->lastValue;
    temp = AdamsMoulton(y_in_a, parameterHigh, a, b, h, temp, 0);
    double highValue = temp->lastValue;

    while ((lowValue - y_in_b) * (highValue - y_in_b) > 0) {
        parameterHigh += 10.0;
        temp = AdamsMoulton(y_in_a, parameterHigh, a, b, h, temp, 0);
        highValue = temp->lastValue;
    }

    double parameter = (parameterLow + parameterHigh) / 2.0;
    while (fabs(y->lastValue - y_in_b) > e) {
        parameter = (parameterLow + parameterHigh) / 2.0;
        temp = AdamsMoulton(y_in_a, parameter, a, b, h, temp, 1);

        if ((temp->lastValue - y_in_b) * (lowValue - y_in_b) < 0) {
            parameterHigh = parameter;
        }
        else {
            parameterLow = parameter;
        }

        y = AdamsMoulton(y_in_a, parameter, a, b, h, y, 1);
        y->n++;
    }
    free(temp);
    return y;
}

Solve* RichardsonExtrapolation(double a, double b, double c, double d, double e, double n) {
    double n1 = n;
    double h1 = (b - a) / (n1 - 1);
    double h2;
    Solve* y1 = NULL;
    Solve* y2 = NULL;
    double error = 0.0;

    do {
        y1 = ShootingMethod(c, d, a, b, h1, e);
        WriteTrueAnswer(a, h1, y1);

        double n2 = 2 * n1;
        h2 = (b - a) / (n2 - 1);
        y2 = ShootingMethod(c, d, a, b, h2, e);
        WriteTrueAnswer(a, h2, y2);

        error = 0.0;
        for (int i = 0; i < y1->length; i++) {
            double diff = fabs(y1->y[i] - y2->y[i * 2]);
            if (diff > error) {
                error = diff;
            }
        }

        if (error > e) {
            free(y1->y);
            free(y1->y_true);
            free(y1);

            n1 *= 2;
            h1 = (b - a) / (n1 - 1);
        }
    } while (error > e);

    return y2;
}

Solve* Solver(double a, double b, double c, double d, double e, double n) {
    Solve* result = RichardsonExtrapolation(a, b, c, d, e, n);
    return result;
}
