/**
 * @ author: Chenyu Bao
 * @ email: bcynuaa@163.com
 * @ date: 2023-11-29 16:41:54
 * @ description:
 */

#ifndef CUBICSPLINE_C_
#define CUBICSPLINE_C_

#include <stdlib.h>

int findFirstIndex(double* x, int n, double x0) {
    int i = 0;
    for (i = 0; i < n-1; i++) {
        if (x0 >= x[i] && x0 <= x[i+1]) {
            break;
        }
    }
    return i;
};

struct CubicSpline {
    int n_; // number of points
    double* x_;
    double* y_;
    double* m_; // second derivative
    double* h_; // interval
    double min_x_;
    double max_x_;
};

double centerDerivative2(double* x, double* y, double* h, int i) {
    return 2 * (-h[i-1]*y[i] + h[i-1]*y[i+1] + h[i]*y[i-1] - h[i]*y[i]) / (h[i-1]*h[i]*(h[i-1]+h[i]));
};

void calInterval(struct CubicSpline* spline) {
    for (int i = 0; i < spline->n_-1; i++) {
        spline->h_[i] = spline->x_[i+1] - spline->x_[i];
    }
    return;
};

void calMatrix(struct CubicSpline* spline, double* a, double* b, double* c, double* d) {
    // regardless the first and last point
    for (int i=0; i < spline->n_-2; i++) {
        a[i] = spline->h_[i];
        b[i] = 2*(spline->h_[i]+spline->h_[i+1]);
        c[i] = spline->h_[i+1];
        d[i] = 6*(spline->y_[i+2]-spline->y_[i+1])/spline->h_[i+1] - 6*(spline->y_[i+1]-spline->y_[i])/spline->h_[i];
    }
    return;
};

void thomas(int n, double* res, double* a, double* b, double* c, double* d) {
    double* l = (double*)malloc(n*sizeof(double));
    double* y = (double*)malloc(n*sizeof(double));
    double* beta = (double*)malloc(n*sizeof(double));
    beta[0] = b[0];
    y[0] = d[0];
    for (int i = 1; i < n; i++) {
        l[i] = a[i]/beta[i-1];
        beta[i] = b[i] - l[i]*c[i-1];
        y[i] = d[i] - l[i]*y[i-1];
    }
    res[n-1] = y[n-1] / beta[n-1];
    for (int i = n-2; i > -1; i--) {
        res[i] = (y[i] - c[i] * res[i+1]) / beta[i];
    }
    free(l); l = NULL;
    free(y); y = NULL;
    free(beta); beta = NULL;
    return;
};

void calCubicSpline(struct CubicSpline* spline) {
    double* a = (double*)malloc( (spline->n_-2)*sizeof(double) );
    double* b = (double*)malloc( (spline->n_-2)*sizeof(double) );
    double* c = (double*)malloc( (spline->n_-2)*sizeof(double) );
    double* d = (double*)malloc( (spline->n_-2)*sizeof(double) );
    double* m = (double*)malloc( (spline->n_-2)*sizeof(double) );
    calMatrix(spline, a, b, c, d);
    d[0] = d[0] - spline->m_[0]*spline->h_[0];
    d[spline->n_-3] = d[spline->n_-3] - spline->m_[spline->n_-1]*spline->h_[spline->n_-2];
    thomas(spline->n_ - 2, m, a, b, c, d);
    for (int i = 0; i < spline->n_-2; i++) {
        spline->m_[i+1] = m[i];
    }
    free(a); a = NULL;
    free(b); b = NULL;
    free(c); c = NULL;
    free(d); d = NULL;
    free(m); m = NULL;
    return;
};

void makeCubicSpline(struct CubicSpline* spline, int n, double* x, double* y) {
    spline->n_ = n;
    spline->x_ = (double*)malloc(n*sizeof(double));
    spline->y_ = (double*)malloc(n*sizeof(double));
    for (int i = 0; i < n; i++) {
        spline->x_[i] = x[i];
        spline->y_[i] = y[i];
    }
    spline->h_ = (double*)malloc((n-1)*sizeof(double));
    spline->m_ = (double*)malloc(n*sizeof(double));
    spline->min_x_ = x[0];
    spline->max_x_ = x[n-1];
    calInterval(spline);
    spline->m_[0] = centerDerivative2(spline->x_, spline->y_, spline->h_, 1);
    spline->m_[n-1] = centerDerivative2(spline->x_, spline->y_, spline->h_, n-2);
    calCubicSpline(spline);
    return;
};

double interpolate(struct CubicSpline* spline, double x) {
    if (x < spline->min_x_) {
        x = spline->min_x_;
        return spline->y_[0] + (spline->y_[1]-spline->y_[0]) * (x-spline->x_[0]) / spline->h_[0];
    }
    else if (x > spline->max_x_) {
        x = spline->max_x_;
        return spline->y_[spline->n_-1] + (spline->y_[spline->n_-1]-spline->y_[spline->n_-2]) * (x-spline->x_[spline->n_-1]) / spline->h_[spline->n_-2];
    }
    else {
        int i = findFirstIndex(spline->x_, spline->n_, x);
        double x0 = spline->x_[i];
        double x1 = spline->x_[i+1];
        double y0 = spline->y_[i];
        double y1 = spline->y_[i+1];
        double m0 = spline->m_[i];
        double m1 = spline->m_[i+1];
        double y = (m0*(x1-x)*(x1-x)*(x1-x) + m1*(x-x0)*(x-x0)*(x-x0)) / (6*spline->h_[i]) + 
            (y0-m0*spline->h_[i]*spline->h_[i]/6)*(x1-x)/spline->h_[i] + 
                (y1-m1*spline->h_[i]*spline->h_[i]/6)*(x-x0)/spline->h_[i];
        return y;
    }
};

double interpolateDerivative1(struct CubicSpline* spline, double x) {
    if (x < spline->min_x_) {
        x = spline->min_x_;
        return (spline->y_[1]-spline->y_[0]) / spline->h_[0];
    }
    else if (x > spline->max_x_) {
        x = spline->max_x_;
        return (spline->y_[spline->n_-1]-spline->y_[spline->n_-2]) / spline->h_[spline->n_-2];
    }
    else {
        int i = findFirstIndex(spline->x_, spline->n_, x);
        double x0 = spline->x_[i];
        double x1 = spline->x_[i+1];
        double y0 = spline->y_[i];
        double y1 = spline->y_[i+1];
        double m0 = spline->m_[i];
        double m1 = spline->m_[i+1];
        double y = (-m0*(x1-x)*(x1-x) + m1*(x-x0)*(x-x0)) / (2*spline->h_[i]) - 
            (y0-m0*spline->h_[i]*spline->h_[i]/6)/spline->h_[i] + 
                (y1-m1*spline->h_[i]*spline->h_[i]/6)/spline->h_[i];
        return y;
    }
};

double interpolateDerivative2(struct CubicSpline* spline, double x) {
    if (x < spline->min_x_) {
        return spline->m_[0];
    }
    else if (x > spline->max_x_) {
        return spline->m_[spline->n_-1];
    }
    else {
        int i = findFirstIndex(spline->x_, spline->n_, x);
        double x0 = spline->x_[i];
        double x1 = spline->x_[i+1];
        double y0 = spline->y_[i];
        double y1 = spline->y_[i+1];
        double m0 = spline->m_[i];
        double m1 = spline->m_[i+1];
        double y = (m0*(x1-x) + m1*(x-x0)) / spline->h_[i];
        return y;
    }
};

void freeCubicSpline(struct CubicSpline* spline) {
    free(spline->x_); spline->x_ = NULL;
    free(spline->y_); spline->y_ = NULL;
    free(spline->m_); spline->m_ = NULL;
    free(spline->h_); spline->h_ = NULL;
    return;
};

#endif // CUBICSPLINE_C_