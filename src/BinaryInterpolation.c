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

/*
  this method is linear interpolation for y-direction
  and cubic spline interpolation for x-direction
  carefully use this method
*/

struct BinaryInterpolation {
    int n_x_;
    int n_y_;
    double* x_;
    double* y_;
    double* z_; // continuous for x-direction,
                // discrete for y-direction
                // z = f(x, y), n_x_ * n_y_
    struct CubicSpline* x_spline_;
    double min_y_;
    double max_y_;
};

void initializeBinaryInterpolation(struct BinaryInterpolation* interpolation, int n_x, int n_y, double* x, double* y, double* z) {
    interpolation->n_x_ = n_x;
    interpolation->n_y_ = n_y;
    interpolation->x_ = x;
    interpolation->y_ = y;
    interpolation->z_ = z;
    interpolation->x_spline_ = (struct CubicSpline*)malloc(n_y*sizeof(struct CubicSpline));
    for (int i = 0; i < n_y; i++) {
        double* z_spline_y = (double*)malloc(n_x*sizeof(double));
        for (int j = 0; j < n_x; j++) {
            z_spline_y[j] = z[i*n_x+j];
        }
        initializeCubicSpline(&interpolation->x_spline_[i], n_x, x, z_spline_y);
        makeCubicSpline(&interpolation->x_spline_[i]);
    }
    return;
};

double binaryInterpolate(struct BinaryInterpolation* interpolation, double x, double y) {
    int i = 0;
    if (y <= interpolation->min_y_) {
        i = 0;
    } else if (y >= interpolation->max_y_) {
        i = interpolation->n_y_ - 2;
    } else {
        for (i = 0; i < interpolation->n_y_-1; i++) {
            if (y >= interpolation->y_[i] && y < interpolation->y_[i+1]) {
                break;
            }
        }
    }
    double z1 = interpolate(&interpolation->x_spline_[i], x);
    double z2 = interpolate(&interpolation->x_spline_[i+1], x);
    return z1 + (z2-z1)/(interpolation->y_[i+1]-interpolation->y_[i])*(y-interpolation->y_[i]);
};

double binaryInterpolateDerivative1(struct BinaryInterpolation* interpolation, double x, double y) {
    int i = 0;
    if (y < interpolation->min_y_) {
        i = 0;
    } else if (y >= interpolation->max_y_) {
        i = interpolation->n_y_ - 2;
    } else {
        for (i = 0; i < interpolation->n_y_-1; i++) {
            if (y >= interpolation->y_[i] && y <= interpolation->y_[i+1]) {
                break;
            }
        }
    }
    double z1 = interpolateDerivative1(&interpolation->x_spline_[i], x);
    double z2 = interpolateDerivative1(&interpolation->x_spline_[i+1], x);
    return z1 + (z2-z1)/(interpolation->y_[i+1]-interpolation->y_[i]);
};

double binaryInterpolateDerivative2(struct BinaryInterpolation* interpolation, double x, double y) {
    int i = 0;
    if (y <= interpolation->min_y_) {
        i = 0;
    } else if (y >= interpolation->max_y_) {
        i = interpolation->n_y_ - 2;
    } else {
        for (i = 0; i < interpolation->n_y_-1; i++) {
            if (y >= interpolation->y_[i] && y <= interpolation->y_[i+1]) {
                break;
            }
        }
    }
    double z1 = interpolateDerivative2(&interpolation->x_spline_[i], x);
    double z2 = interpolateDerivative2(&interpolation->x_spline_[i+1], x);
    return z1 + (z2-z1)/(interpolation->y_[i+1]-interpolation->y_[i]);
};

#endif // BINARYINTERPOLATION_C_