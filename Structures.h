#ifndef STRUCT_H
#define STRUCT_H

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
	double operator*(const Point& other) const { return x * other.x + y * other.y; }

};

struct Matrix {
	Point row_1;
	Point row_2;

	double det() const
	{
		return row_1.x * row_2.y - row_1.y * row_2.x;
	}

	Matrix inverse() const {
		double determinant = det();

		double inv_det = 1.0 / determinant;
		return Matrix{
			Point{ row_2.y, -row_1.y } * inv_det,
			Point{ -row_2.x, row_1.x } * inv_det
		};
	}

	Point operator*(const Point& p) const {
		return {
			row_1.x * p.x + row_1.y * p.y,
			row_2.x * p.x + row_2.y * p.y
		};
	}
};

#endif
