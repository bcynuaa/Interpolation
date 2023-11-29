/**
 * @ author: Chenyu Bao
 * @ email: bcynuaa@163.com
 * @ date: 2023-11-29 20:30:40
 * @ description:
 */

#include <stdio.h>
#include "CubicSpline.c"

int main() {
    int n = 11;
    double x[11] = {0. , 0.5, 1. , 1.5, 2. , 2.5, 3. , 3.5, 4. , 4.5, 5. };
    double y[11] = {-1.        , -0.60651042, -0.36666667, -0.21015625, -0.06666667,
        0.16536458,  0.65      ,  1.64505208,  3.53333333,  6.85390625, 12.33333333};
    struct CubicSpline spline;
    makeCubicSpline(&spline, n, x, y);
    int i = 0;
    int minor_n = 101;
    double* x_minor = (double*)malloc(minor_n*sizeof(double));
    for (i = 0; i < minor_n; i++) {
        x_minor[i] = 5 / (double)(minor_n - 1) * i;
    }
    double* y_minor = (double*)malloc(minor_n*sizeof(double));
    for (i = 0; i < minor_n; i++) {
        y_minor[i] = interpolate(&spline, x_minor[i]);
    }
    double* y1_minor = (double*)malloc(minor_n*sizeof(double));
    for (i = 0; i < minor_n; i++) {
        y1_minor[i] = interpolateDerivative1(&spline, x_minor[i]);
    }
    double* y2_minor = (double*)malloc(minor_n*sizeof(double));
    for (i = 0; i < minor_n; i++) {
        y2_minor[i] = interpolateDerivative2(&spline, x_minor[i]);
    }
    // write to file
    FILE* fp;
    fp = fopen("CubicSplineTest.dat", "w");
    for (i = 0; i < minor_n; i++) {
        fprintf(fp, "%f %f %f %f\n", x_minor[i], y_minor[i], y1_minor[i], y2_minor[i]);
    }
    fclose(fp);
    freeCubicSpline(&spline);
    return 0;
}