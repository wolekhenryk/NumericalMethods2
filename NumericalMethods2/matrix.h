#pragma once
#include <iostream>
#include <cmath>
#include <vector>

class matrix {
	int rows_, cols_;
	std::vector<std::vector<double>> data_;

public:
	matrix(int rows, int cols);
	explicit matrix(const std::vector<std::vector<double>>& data);

	[[nodiscard]] int rows() const;
	[[nodiscard]] int cols() const;

	// Access element
	double& operator()(int row, int col);

	const double& operator()(int row, int col) const;

	[[nodiscard]] matrix operator*(const matrix& other) const;

	[[nodiscard]] matrix inverse() const;

	void fill_five_diagonals(double first, double second, double third);

	// Print matrix
	friend std::ostream& operator<<(std::ostream& os, const matrix& m);
};
