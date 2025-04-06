#ifndef METHOD_H
#define METHOD_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>

using namespace std;

struct Point {
    double x, y;

    Point operator+(const Point& other) const { return { x + other.x, y + other.y }; }
    Point operator-(const Point& other) const { return { x - other.x, y - other.y }; }
    Point operator*(double scalar) const { return { x * scalar, y * scalar }; }

    // Methods
    double norm() const { return sqrt(x * x + y * y); }
};


double f(const Point& p) {
    return 0.5 * p.x * p.x + 0.5 * p.y * p.y + sin(p.x) * sin(p.x) * cos(p.y) * cos(p.y);
}


Point gradient(const Point& p) {
    double df_dx = p.x + 2.0 * sin(p.x) * cos(p.x) * cos(p.y) * cos(p.y);
    double df_dy = p.y - 2.0 * sin(p.x) * sin(p.x) * cos(p.y) * sin(p.y);
    return { df_dx, df_dy };
}

// One-dimensional minimization method for finding the optimal step
double golden_section_search(const Point& x, const Point& d, double a, double b, double eps, int max_iter = 1000) {
    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    double c1 = b - (b - a) / phi;
    double ñ2 = a + (b - a) / phi;
    int it = 0;

    while (fabs(c1 - ñ2) > eps) {
        if (it == max_iter) {
            cout << "The maximum number of iterations has been reached in an inner method." << endl;
            break;
        }

        if (f(x + d * c1) < f(x + d * ñ2)) {
            b = ñ2;
        }
        else {
            a = c1;
        }
        c1 = b - (b - a) / phi;
        ñ2 = a + (b - a) / phi;

        it++;
    }

    return (a + b) / 2.0;
}

Point steepest_descent(Point start_point, double L, double eps, int max_iter = 1000) {
    Point point = start_point;
    Point prev_point; // To check orthogonality
    int k = -1; // Thtat's kinda bad...

    while (gradient(point).norm() > eps) {

        k++;

        if (k == max_iter) {
            cout << "The maximum number of iterations has been reached in an gradient method." << endl;
            break;
        }

        Point grad = gradient(point);
        if (grad.norm() < eps) {
            break;
        }

        // Antigradient direction
        Point d = grad * (-1.0);

        // Step optimization using the golden section method
        double alpha = golden_section_search(point, d, 0.0, (double)2 / L, eps);

        // Update point
        point = point + d * alpha;

        // Checking the orthogonality of gradients at adjacent iterations
        if (k >= 1) {
            Point prev_grad = gradient(prev_point);
            double dot_product = prev_grad.x * grad.x + prev_grad.y * grad.y;
            cout << "Step " << k << ": Dot product of gradients = " << dot_product << endl;
        }
        prev_point = point;
    }

    cout << "\nNumber of iterations of steepest descent method: " << k << endl;

    return point;
}

#endif