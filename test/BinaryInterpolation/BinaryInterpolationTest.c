/**
 * @ author: Chenyu Bao
 * @ email: bcynuaa@163.com
 * @ date: 2023-11-29 21:23:20
 * @ description:
 */

#include <stdio.h>
#include "BinaryInterpolation.c"

int main() {
    struct BinaryInterpolation bin_ip;
    int nx = 11;
    int ny = 6;
    double x[11] = {0. , 0.5, 1. , 1.5, 2. , 2.5, 3. , 3.5, 4. , 4.5, 5. };
    double y[6] = {-1.,  0.,  1.,  2.,  3.,  4.};
    double z[121] = {2.00000000e+00,  1.21302083e+00,  7.33333333e-01,  4.20312500e-01,
        1.33333333e-01, -3.30729167e-01, -1.30000000e+00, -3.29010417e+00,
       -7.06666667e+00, -1.37078125e+01, -2.46666667e+01, -1.00000000e+00,
       -6.06510417e-01, -3.66666667e-01, -2.10156250e-01, -6.66666667e-02,
        1.65364583e-01,  6.50000000e-01,  1.64505208e+00,  3.53333333e+00,
        6.85390625e+00,  1.23333333e+01, -4.00000000e+00, -2.42604167e+00,
       -1.46666667e+00, -8.40625000e-01, -2.66666667e-01,  6.61458333e-01,
        2.60000000e+00,  6.58020833e+00,  1.41333333e+01,  2.74156250e+01,
        4.93333333e+01, -7.00000000e+00, -4.24557292e+00, -2.56666667e+00,
       -1.47109375e+00, -4.66666667e-01,  1.15755208e+00,  4.55000000e+00,
        1.15153646e+01,  2.47333333e+01,  4.79773437e+01,  8.63333333e+01,
       -1.00000000e+01, -6.06510417e+00, -3.66666667e+00, -2.10156250e+00,
       -6.66666667e-01,  1.65364583e+00,  6.50000000e+00,  1.64505208e+01,
        3.53333333e+01,  6.85390625e+01,  1.23333333e+02, -1.30000000e+01,
       -7.88463542e+00, -4.76666667e+00, -2.73203125e+00, -8.66666667e-01,
        2.14973958e+00,  8.45000000e+00,  2.13856771e+01,  4.59333333e+01,
        8.91007812e+01,  1.60333333e+02};
    makeBinaryInterpolation(&bin_ip, nx, ny, x, y, z);
    int minor_n = 101;
    double* minor_x = (double*)malloc(minor_n*sizeof(double));
    double* minor_y = (double*)malloc(minor_n*sizeof(double));
    double* minor_z = (double*)malloc(minor_n*minor_n*sizeof(double));
    double x_begin = 0.;
    double x_end = 5.;
    double y_begin = -1.;
    double y_end = 4.;
    for (int i = 0; i < minor_n; i++) {
        minor_x[i] = x_begin + (x_end - x_begin) / (minor_n - 1) * i;
        minor_y[i] = y_begin + (y_end - y_begin) / (minor_n - 1) * i;
    }
    for (int i = 0; i < minor_n; i++) {
        for (int j = 0; j < minor_n; j++) {
            minor_z[i*minor_n+j] = binaryInterpolate(&bin_ip, minor_x[j], minor_y[i]);
        }
    }
    FILE* fp = fopen("BinaryInterpolationTest.dat", "w");
    for (int i = 0; i < minor_n; i++) {
        for (int j = 0; j < minor_n; j++) {
            fprintf(fp, "%lf ", minor_z[i*minor_n+j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}