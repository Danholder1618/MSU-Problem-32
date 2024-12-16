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
    return 1 - x * x + 2 * x - 11 * x + (x * x * x * x) / 6.0 + (x * x) / 2.0;
}

Solve* AdamsMoulton6RangSolver(double y1_0, double y2_0, double a, double b, double h, Solve* solve, int needSolveWithPoints)
{
    unsigned long long pointsCount = fabs(b - a) / h + 1;
    double* y1 = (double*)malloc(pointsCount * sizeof(double));
    double* y2 = (double*)malloc(pointsCount * sizeof(double));

    y1[0] = y1_0;
    y2[0] = y2_0;

    double t = a;

    for (int i = 1; i < 6; i++) {
        double k1_y1 = h * F1(y2[i - 1]);
        double k1_y2 = h * F2(t, y1[i - 1], y2[i - 1]);

        double k2_y1 = h * F1(y2[i - 1] + k1_y2 / 2);
        double k2_y2 = h * F2(t + h / 2, y1[i - 1] + k1_y1 / 2, y2[i - 1] + k1_y2 / 2);

        double k3_y1 = h * F1(y2[i - 1] + k2_y2 / 2);
        double k3_y2 = h * F2(t + h / 2, y1[i - 1] + k2_y1 / 2, y2[i - 1] + k2_y2 / 2);

        double k4_y1 = h * F1(y2[i - 1] + k3_y2);
        double k4_y2 = h * F2(t + h, y1[i - 1] + k3_y1, y2[i - 1] + k3_y2);

        y1[i] = y1[i - 1] + (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1) / 6;
        y2[i] = y2[i - 1] + (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2) / 6;

        t += h;
    }

    for (int i = 5; i < pointsCount - 1; i++) {
        double predictor_y1 = y1[i] + h * F1(y2[i]);
        double predictor_y2 = y2[i] + h * F2(t, y1[i], y2[i]);

        double corrector_y1, corrector_y2;
        for (int iter = 0; iter < 3; iter++) {
            corrector_y1 = y1[i] + h / 720.0 * (1901 * F1(y2[i]) - 2774 * F1(y2[i - 1]) + 2616 * F1(y2[i - 2]) - 1274 * F1(y2[i - 3]) + 251 * F1(y2[i - 4]));
            corrector_y2 = y2[i] + h / 720.0 * (1901 * F2(t, y1[i], y2[i]) - 2774 * F2(t - h, y1[i - 1], y2[i - 1]) + 2616 * F2(t - 2 * h, y1[i - 2], y2[i - 2]) - 1274 * F2(t - 3 * h, y1[i - 3], y2[i - 3]) + 251 * F2(t - 4 * h, y1[i - 4], y2[i - 4]));

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

Solve* ShootingMethod(double y_in_a, double y_in_b, double a, double b, double h, double e,
    Solve* AdamsMoulton6RangSolver(double y1_0, double y2_0, double a, double b, double h,
        Solve* solve, int needSolveWithPoints))
{
    Solve* y = (Solve*)malloc(sizeof(Solve));
    y->n = 0;
    Solve* temp = (Solve*)malloc(sizeof(Solve));

    double parameterLow = -10.0;
    double parameterHigh = 10.0;

    temp = AdamsMoulton6RangSolver(y_in_a, parameterLow, a, b, h, temp, 0);
    double lowValue = temp->lastValue;
    temp = AdamsMoulton6RangSolver(y_in_a, parameterHigh, a, b, h, temp, 0);
    double highValue = temp->lastValue;

    while ((lowValue - y_in_b) * (highValue - y_in_b) > 0) {
        parameterHigh += 10.0;
        temp = AdamsMoulton6RangSolver(y_in_a, parameterHigh, a, b, h, temp, 0);
        highValue = temp->lastValue;
    }

    double parameter = (parameterLow + parameterHigh) / 2.0;
    while (fabs(y->lastValue - y_in_b) > e) {
        parameter = (parameterLow + parameterHigh) / 2.0;
        temp = AdamsMoulton6RangSolver(y_in_a, parameter, a, b, h, temp, 1);

        if ((temp->lastValue - y_in_b) * (lowValue - y_in_b) < 0) {
            parameterHigh = parameter;
        }
        else {
            parameterLow = parameter;
        }

        y = AdamsMoulton6RangSolver(y_in_a, parameter, a, b, h, y, 1);
        y->n++;
    }
    free(temp);
    return y;
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

Solve* Solver(double a, double b, double c, double d, double h, double e)
{
    Solve* y = ShootingMethod(c, d, a, b, h, e, AdamsMoulton6RangSolver);
    WriteTrueAnswer(a, h, y);
    return y;
}