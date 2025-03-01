#include "pbPlots.h"
#include "supportLib.h"
#include "header.h"

#define _CRT_SECURE_NO_WARNINGS

typedef struct InputData {
    double a;
    double b;
    double c;
    double d;
    double e;
    double n;
} InputData;

InputData* ReadInput() {
    FILE* input = fopen("input.txt", "r");
    if (input == NULL) {
        printf("Error opening file\n");
        return NULL;
    }

    InputData* inputData = (InputData*)malloc(sizeof(InputData));
    int values = 0;
    values += fscanf(input, "%lf", &inputData->a);
    values += fscanf(input, "%lf", &inputData->b);
    values += fscanf(input, "%lf", &inputData->c);
    values += fscanf(input, "%lf", &inputData->d);
    values += fscanf(input, "%lf", &inputData->e);
    values += fscanf(input, "%lf", &inputData->n);

    fclose(input);

    if (values != 6) {
        printf("Invalid number of values in input file. Expected 6.\n");
        free(inputData);
        return NULL;
    }

    return inputData;
}

_Bool WriteOutputData(Solve* solve) {
    FILE* output = fopen("output.txt", "w");
    if (output == NULL) return false;

    fprintf(output, "%d\n", solve->n);
    for (int i = 0; i < solve->length; i++) {
        fprintf(output, "%.10lf ", solve->y[i]);
    }
    fclose(output);
    return true;
}

_Bool DrawGraphics(double a, double h, Solve* solve) {
    double* xs = (double*)malloc(solve->length * sizeof(double));
    double x = a;
    double width = 1280;
    double high = 720;

    for (int i = 0; i < solve->length; i++) {
        xs[i] = x;
        x += h;
    }

    double* y = solve->y;
    double* y_true = solve->y_true;

    _Bool success;
    ScatterPlotSeries* y_plot = GetDefaultScatterPlotSeriesSettings();
    y_plot->linearInterpolation = true;
    y_plot->xs = xs;
    y_plot->xsLength = solve->length;
    y_plot->ys = y;
    y_plot->ysLength = solve->length;
    y_plot->lineType = L"solid";
    y_plot->lineTypeLength = wcslen(y_plot->lineType);
    y_plot->lineThickness = 2;
    y_plot->color = CreateRGBColor(0.7, 0.05, 0.05);

    ScatterPlotSeries* y_true_plot = GetDefaultScatterPlotSeriesSettings();
    y_true_plot->linearInterpolation = true;
    y_true_plot->xs = xs;
    y_true_plot->xsLength = solve->length;
    y_true_plot->ys = y_true;
    y_true_plot->ysLength = solve->length;
    y_true_plot->lineType = L"solid";
    y_true_plot->lineTypeLength = wcslen(y_true_plot->lineType);
    y_true_plot->lineThickness = 2;
    y_true_plot->color = CreateRGBColor(0.05, 0.7, 0.05);

    ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();
    settings->width = width;
    settings->height = high;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = L"Graphic of true solution and its approximation";
    settings->titleLength = wcslen(settings->title);
    settings->xLabel = L"X axis";
    settings->xLabelLength = wcslen(settings->xLabel);
    settings->yLabel = L"Y axis";
    settings->yLabelLength = wcslen(settings->yLabel);

    ScatterPlotSeries* s[] = { y_plot, y_true_plot };
    settings->scatterPlotSeries = s;
    settings->scatterPlotSeriesLength = 2;

    RGBABitmapImageReference* canvasReference = CreateRGBABitmapImageReference();
    StringReference* errorMessage = (StringReference*)malloc(sizeof(StringReference));
    success = DrawScatterPlotFromSettings(canvasReference, settings, errorMessage);

    if (success) {
        size_t length;
        double* pngdata = ConvertToPNG(&length, canvasReference->image);
        WriteToFile(pngdata, length, "SolveGraphic.png");
        DeleteImage(canvasReference->image);
    }
    else {
        fprintf(stderr, "Error: ");
        for (int i = 0; i < errorMessage->stringLength; i++) {
            fprintf(stderr, "%c", errorMessage->string[i]);
        }
        fprintf(stderr, "\n");
    }

    return success ? 0 : 1;
}

int main() {
    InputData* inputData = ReadInput();
    if (inputData == NULL) {
        return 1;
    }

    Solve* solve = Solver(inputData->a, inputData->b, inputData->c, inputData->d, inputData->e, inputData->n);

    double h = (inputData->b - inputData->a) / (solve->length - 1);
    DrawGraphics(inputData->a, h, solve);

    if (WriteOutputData(solve)) {
        printf("Process Done\n");
    }

    free(inputData);
    free(solve->y);
    free(solve->y_true);
    free(solve);

    return 0;
}
