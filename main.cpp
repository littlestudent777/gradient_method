#include "Method.h"

int main() {
    Point opt_point = { 0.0, 0.0 }; // Actual min point
    double opt_val = 0.0; // Actual min value
    Point start_point = { 1.0, 1.0 }; // Initial point
    double L = 4;

    Point res_point = steepest_descent(start_point, L, 1e-6);

    cout << "\nOptimal point: (" << res_point.x << ", " << res_point.y << ")" << endl;
    cout << "Optimal value: " << f(res_point) << endl;

    return 0;
}
