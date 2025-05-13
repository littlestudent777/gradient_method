#ifndef METHOD_H
#define METHOD_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>
#include "Structures.h"

using namespace std;

double f(const Point& p) {
    return 0.5 * p.x * p.x + 0.5 * p.y * p.y + sin(p.x) * sin(p.x) * cos(p.y) * cos(p.y);
}

Point gradient(const Point& p) {
    double df_dx = p.x + 2.0 * sin(p.x) * cos(p.x) * cos(p.y) * cos(p.y);
    double df_dy = p.y - 2.0 * sin(p.x) * sin(p.x) * cos(p.y) * sin(p.y);
    return { df_dx, df_dy };
}

Matrix hessian(const Point& p) {
    Matrix res;
    double cos2x = cos(2 * p.x);
    double cos2y = cos(2 * p.y);
    double sin2x = sin(2 * p.x);
    double sin2y = sin(2 * p.y);
    double cosy2 = pow(cos(p.y), 2);
    double sinx2 = pow(sin(p.x), 2);

    res.row_1 = { 1 + 2 * cos2x * cosy2, -sin2x * sin2y };
    res.row_2 = { -sin2x * sin2y,      1 - 2 * sinx2 * cos2y };

    return res;
}

void find_interval(const Point& x, const Point& d, double& a, double& b, double alpha0 = 0.5, double gamma = 2.0) {
    double alpha_prev = 0.0;
    double alpha_curr = alpha0;
    double phi_prev = f(x + d * alpha_prev);
    double phi_curr = f(x + d * alpha_curr);

    // If the function immediately increases, we take [0, alpha0]
    if (phi_curr > phi_prev) {
        a = 0.0;
        b = alpha_curr;
        return;
    }

    // Expand the interval until the function starts to increase
    while (phi_curr < phi_prev) {
        alpha_prev = alpha_curr;
        phi_prev = phi_curr;
        alpha_curr *= gamma;
        phi_curr = f(x + d * alpha_curr);
    }

    a = alpha_prev / gamma;
    b = alpha_curr;
}

// One-dimensional minimization method for finding the optimal step
double golden_section_search(const Point& x, const Point& d, double a, double b, double eps, int max_iter = 1000) {
    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    double c1 = b - (b - a) / phi;
    double c2 = a + (b - a) / phi;
    int it = 0;

    while (fabs(c1 - c2) > eps) {
        if (it == max_iter) {
            cout << "The maximum number of iterations has been reached in an inner method." << endl;
            break;
        }

        if (f(x + d * c1) < f(x + d * c2)) {
            b = c2;
        }
        else {
            a = c1;
        }
        c1 = b - (b - a) / phi;
        c2 = a + (b - a) / phi;

        it++;
    }

    return (a + b) / 2.0;
}


Point steepest_descent(Point start_point, double eps, int max_iter = 1000) {
    Point point = start_point;
    Point prev_point; // To check orthogonality
    int it = 0;

    for (int k = 0; k < max_iter; ++k) {

        Point grad = gradient(point);

        if (grad.norm() < eps) {
            it = k;
            break;
        }

        // Antigradient direction
         Point d = grad * (-1.0);

         double a, b;
         find_interval(point, d, a, b);

        // Step optimization using the golden section method
        double alpha = golden_section_search(point, d, a, b, eps);

        // Update point
        point = point + d * alpha;

        //cout << "k = " << k << ", ||gradf(x_k)||/m = " << gradient(point).norm()/0.55 << endl;
        //cout << "k = " << k << ", ||x_k - x*|| = " << point.norm() << endl;
        //cout << "k = " << k << ", x_k : (" << point.x << ", " << point.y << ")" << endl;

        // Checking the orthogonality of gradients at adjacent iterations
        if (k == 0) {
            cout << "Dot product of gradients in steepest descent method:" << endl;
        }
        else {
            Point prev_grad = gradient(prev_point);
            double dot_product = prev_grad.x * grad.x + prev_grad.y * grad.y;
            cout << "Step " << k << ": " << dot_product << endl;
        }
        prev_point = point;
    }

    cout << "\nNumber of iterations of steepest descent method: " << it << endl;

    return point;
}

Point compute_D(const Point& p) {
    return {
        1.0 / (1.0 + 2.0 * cos(p.y) * cos(p.y)),  // D11
        1.0 / (1.0 + 2.0 * sin(p.x) * sin(p.x))   // D22
    };
}

// Custom method
Point bottom_search_method(const Point& x0, double eps, int max_iter = 1000) {
    Point x = x0;
    int it = 0;
    for (int k = 0; k < max_iter; ++k) {
        Point g = gradient(x);

        if (g.norm() < eps) {
            it = k;
            break;
        }

        Point D = compute_D(x);

        // Direction: d = -D * g (element-wise multiplication)
        Point d = { -D.x * g.x, -D.y * g.y };

        double a, b;
        //find_interval(x, d, a, b);

        // Step optimization using the golden section method
        double alpha = golden_section_search(x, d, 0.0, 1.0, eps);

        x = x + d * alpha;

        // cout << "k = " << k << ", x_k : (" << x.x << ", " << x.y << ")" << endl;

    }
    
    cout << "\nNumber of iterations on custom: " << it << endl;
    return x;
}

Point pschenychny(const Point& x0, double epsilon, int maxIterations = 1000, double c = 0.1, double tau = 0.5) {
    Point x = x0;
    int it = 0;

    for (int i = 0; i < maxIterations; ++i) {
        Point grad = gradient(x);
        if (grad.norm() < epsilon) {
            it = i;
            break;
        }

        Matrix H = hessian(x);
        Point direction;

        // If the det of Hessian matrix is close to zero, use gradient descent
        if (fabs(H.det()) < 1e-10) {
            direction = grad * (-1.0);
        }
        else {
            H = H.inverse();
            direction = H * grad * (-1.0);
        }

        // Check that the direction is the direction of descent
        double grad_dot_dir = grad * direction;
        if (grad_dot_dir >= 0) {
            direction = grad * (-1.0);
            grad_dot_dir = grad * direction;
        }

        // Step selection
        double alpha = 1.0;
        double f_cur = f(x);
        bool IsStepFound = false;

        for (int step = 0; step < 100; ++step) {
            Point x_new = x + direction * alpha;
            double f_new = f(x_new);

            // Checking the Armijo condition (sufficient decrease)
            if (f_new <= f_cur + c * alpha * grad_dot_dir) {
                x = x_new;
                IsStepFound = true;
                break;
            }
            // If the condition is not met, do step splitting
            alpha *= tau;
        }

        // If we haven't found a suitable step after 100 attempts, 
        // we take the step with the minimum alpha found
        if (!IsStepFound) {
            x = x + direction * alpha;
        }

        // cout << "k = " << i << ", x_k : (" << x.x << ", " << x.y << ")" << endl;

    }

    cout << "\nNumber of iterations on second order method: " << it << endl;
    return x;
}

#endif