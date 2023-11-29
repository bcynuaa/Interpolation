/**
 * @ author: Chenyu Bao
 * @ email: bcynuaa@163.com
 * @ date: 2023-11-29 21:07:51
 * @ description:
 */

#ifndef BINARYINTERPOLATION_C_
#define BINARYINTERPOLATION_C_

#include <stdlib.h>
#include "CubicSpline.c"

struct BinaryInterpolation {
    int n_x_;
    int n_y_;
    double* x_;
    double* y_;
    struct CubicSpline* x_spline_;
    double min_y_;
    double max_y_;
};

void makeBinaryInterpolation(struct BinaryInterpolation* interpolation, int n_x, int n_y, double* x, double* y, double* z) {
    interpolation->n_x_ = n_x;
    interpolation->n_y_ = n_y;
    interpolation->x_ = (double*)malloc(n_x*sizeof(double));
    interpolation->y_ = (double*)malloc(n_y*sizeof(double));
    for (int i = 0; i < n_x; i++) {
        interpolation->x_[i] = x[i];
    }
    for (int i = 0; i < n_y; i++) {
        interpolation->y_[i] = y[i];
    }
    interpolation->min_y_ = y[0];
    interpolation->max_y_ = y[n_y-1];
    interpolation->x_spline_ = (struct CubicSpline*)malloc(n_y*sizeof(struct CubicSpline));
    double* z_spline_y = (double*)malloc(n_x*sizeof(double));
    for (int i = 0; i < n_y; i++) {
        for (int j = 0; j < n_x; j++) {
            z_spline_y[j] = z[i*n_x+j];
        }
        makeCubicSpline(&interpolation->x_spline_[i], n_x, x, z_spline_y);
    }
    free(z_spline_y); z_spline_y = NULL;
    return;
};

double binaryInterpolate(struct BinaryInterpolation* interpolation, double x, double y) {
    double* z_spline_y = (double*)malloc(interpolation->n_y_*sizeof(double));
    for (int i = 0; i < interpolation->n_y_; i++) {
        z_spline_y[i] = interpolate(&interpolation->x_spline_[i], x);
    }
    struct CubicSpline y_spline;
    makeCubicSpline(&y_spline, interpolation->n_y_, interpolation->y_, z_spline_y);
    double res = interpolate(&y_spline, y);
    free(z_spline_y); z_spline_y = NULL;
    freeCubicSpline(&y_spline);
    return res;
};

double binaryInterpolateDerivativeX1(struct BinaryInterpolation* interpolation, double x, double y) {
    double* z_spline_y = (double*)malloc(interpolation->n_y_*sizeof(double));
    for (int i = 0; i < interpolation->n_y_; i++) {
        z_spline_y[i] = interpolateDerivative1(&interpolation->x_spline_[i], x);
    }
    struct CubicSpline y_spline;
    makeCubicSpline(&y_spline, interpolation->n_y_, interpolation->y_, z_spline_y);
    double res = interpolate(&y_spline, y);
    free(z_spline_y); z_spline_y = NULL;
    freeCubicSpline(&y_spline);
    return res;
};

double binaryInterpolateDerivativeX2(struct BinaryInterpolation* interpolation, double x, double y) {
    double* z_spline_y = (double*)malloc(interpolation->n_y_*sizeof(double));
    for (int i = 0; i < interpolation->n_y_; i++) {
        z_spline_y[i] = interpolateDerivative2(&interpolation->x_spline_[i], x);
    }
    struct CubicSpline y_spline;
    makeCubicSpline(&y_spline, interpolation->n_y_, interpolation->y_, z_spline_y);
    double res = interpolate(&y_spline, y);
    free(z_spline_y); z_spline_y = NULL;
    freeCubicSpline(&y_spline);
    return res;
};

double binaryInterpolateDerivativeY1(struct BinaryInterpolation* interpolation, double x, double y) {
    double* z_spline_y = (double*)malloc(interpolation->n_y_*sizeof(double));
    for (int i = 0; i < interpolation->n_y_; i++) {
        z_spline_y[i] = interpolate(&interpolation->x_spline_[i], x);
    }
    struct CubicSpline y_spline;
    makeCubicSpline(&y_spline, interpolation->n_y_, interpolation->y_, z_spline_y);
    double res = interpolateDerivative1(&y_spline, y);
    free(z_spline_y); z_spline_y = NULL;
    freeCubicSpline(&y_spline);
    return res;
};

double binaryInterpolateDerivativeY2(struct BinaryInterpolation* interpolation, double x, double y) {
    double* z_spline_y = (double*)malloc(interpolation->n_y_*sizeof(double));
    for (int i = 0; i < interpolation->n_y_; i++) {
        z_spline_y[i] = interpolate(&interpolation->x_spline_[i], x);
    }
    struct CubicSpline y_spline;
    makeCubicSpline(&y_spline, interpolation->n_y_, interpolation->y_, z_spline_y);
    double res = interpolateDerivative2(&y_spline, y);
    free(z_spline_y); z_spline_y = NULL;
    freeCubicSpline(&y_spline);
    return res;
};

/*
  below is the old version,
  linear interpolation is used in y-direction
  and cubic spline interpolation is used in x-direction
  to speed up the interpolation, you may need the old version
  to precise the interpolation, you may need the new version
*/

// double binaryInterpolate(struct BinaryInterpolation* interpolation, double x, double y) {
//     int i = 0;
//     if (y <= interpolation->min_y_) {
//         i = 0;
//     } else if (y >= interpolation->max_y_) {
//         i = interpolation->n_y_ - 2;
//     } else {
//         for (i = 0; i < interpolation->n_y_-1; i++) {
//             if (y >= interpolation->y_[i] && y < interpolation->y_[i+1]) {
//                 break;
//             }
//         }
//     }
//     double z1 = interpolate(&interpolation->x_spline_[i], x);
//     double z2 = interpolate(&interpolation->x_spline_[i+1], x);
//     return z1 + (z2-z1)/(interpolation->y_[i+1]-interpolation->y_[i])*(y-interpolation->y_[i]);
// };

// double binaryInterpolateDerivative1(struct BinaryInterpolation* interpolation, double x, double y) {
//     int i = 0;
//     if (y < interpolation->min_y_) {
//         i = 0;
//     } else if (y >= interpolation->max_y_) {
//         i = interpolation->n_y_ - 2;
//     } else {
//         for (i = 0; i < interpolation->n_y_-1; i++) {
//             if (y >= interpolation->y_[i] && y <= interpolation->y_[i+1]) {
//                 break;
//             }
//         }
//     }
//     double z1 = interpolateDerivative1(&interpolation->x_spline_[i], x);
//     double z2 = interpolateDerivative1(&interpolation->x_spline_[i+1], x);
//     return z1 + (z2-z1)/(interpolation->y_[i+1]-interpolation->y_[i]);
// };

// double binaryInterpolateDerivative2(struct BinaryInterpolation* interpolation, double x, double y) {
//     int i = 0;
//     if (y <= interpolation->min_y_) {
//         i = 0;
//     } else if (y >= interpolation->max_y_) {
//         i = interpolation->n_y_ - 2;
//     } else {
//         for (i = 0; i < interpolation->n_y_-1; i++) {
//             if (y >= interpolation->y_[i] && y <= interpolation->y_[i+1]) {
//                 break;
//             }
//         }
//     }
//     double z1 = interpolateDerivative2(&interpolation->x_spline_[i], x);
//     double z2 = interpolateDerivative2(&interpolation->x_spline_[i+1], x);
//     return z1 + (z2-z1)/(interpolation->y_[i+1]-interpolation->y_[i]);
// };

#endif // BINARYINTERPOLATION_C_