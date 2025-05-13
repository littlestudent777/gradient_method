#include "Method.h"

// Counting error according to the first norm
double err(const Point& a, const Point& b) {
    return fabs(a.x - b.x) + fabs(a.y - b.y);
}

int main() {
    Point opt_point = { 0.0, 0.0 }; // Actual min point
    double opt_val = 0.0; // Actual min value
    Point start_point = { 1.0, 1.0 }; // Initial point

    double eps = 1e-5;

    cout << "Steepest descent: " << endl;

    Point res_point = steepest_descent(start_point, eps);
    double res_value = f(res_point);

    cout << "Optimal point: (" << res_point.x << ", " << res_point.y << ")" << endl;
    cout << "Min value: " << res_value << endl << endl;


    cout << "Custom method: " << endl;

    Point res_p_new = bottom_search_method(start_point, eps);
    double res_val_new = f(res_p_new);

    cout << "Optimal point: x = (" << res_p_new.x << ", " << res_p_new.y << ")" << endl;
    cout << "Function value: f(x) = " << res_val_new << endl << endl;

   // cout << "\nOptimal point error: " << err(opt_point, x_opt) << endl;
   // cout << "\nMin value error: " << fabs(opt_val - res_val_new) << endl;

    cout << "Pschenychny method: " << endl;

    Point sec_or_res = pschenychny(start_point, eps);
    cout << "\nOptimal point: x = (" << sec_or_res.x << ", " << sec_or_res.y << ")" << endl;
    cout << "Function value: f(x) = " << f(sec_or_res) << endl;

    return 0;
}
