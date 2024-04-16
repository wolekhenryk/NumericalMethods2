#include "matrix.h"

#include <algorithm>
#include <stdexcept>

matrix::matrix(const int rows, const int cols): rows_(rows), cols_(cols), data_(rows, std::vector<double>(cols, 0)) {
}

matrix::matrix(const std::vector<std::vector<double>>& data) {
	rows_ = static_cast<int>(data.size());
	cols_ = static_cast<int>(data[0].size());
	data_ = data;
}

int matrix::rows() const { return rows_; }

int matrix::cols() const { return cols_; }

double& matrix::operator()(const int row, const int col) {
	return data_[row][col];
}

const double& matrix::operator()(const int row, const int col) const {
	return data_[row][col];
}

matrix matrix::operator*(const matrix& other) const {
	if (cols_ != other.rows_) throw std::invalid_argument("Matrix dimensions must agree.");
	matrix result(rows_, other.cols_);

	#pragma omp parallel for
	for (int i = 0; i < rows_; ++i) {
		for (int k = 0; k < cols_; ++k) { // Innermost loop changes from j to k
			const auto temp = (*this)(i, k);
			for (int j = 0; j < other.cols_; ++j) {
				result(i, j) += temp * other(k, j);
			}
		}
	}
	return result;
}

matrix matrix::inverse() const {
	if (rows_ != cols_) throw std::invalid_argument("Only square matrices can be inverted.");

	int n = rows_;
	matrix result(n, n * 2);
	// Create the augmented matrix [A|I]
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result(i, j) = (*this)(i, j);
		}
		result(i, n + i) = 1; // Place identity matrix on the right
	}

	// Perform Gaussian elimination to convert [A|I] to [I|A^-1]
	for (int i = 0; i < n; ++i) {
		// Make the diagonal element 1
		double diagonal = result(i, i);
		if (std::fabs(diagonal) < 1e-10) {
			// checking for zero
			// Need to swap with a non-zero row
			int swap_row = i + 1;
			while (swap_row < n && std::fabs(result(swap_row, i)) < 1e-10)
				++swap_row;
			if (swap_row == n)
				throw std::runtime_error("Matrix is singular and cannot be inverted.");
			// Swap the rows
			std::swap_ranges(result.data_[i].begin(), result.data_[i].end(), result.data_[swap_row].begin());
			diagonal = result(i, i);
		}
		for (int j = 0; j < n * 2; ++j) {
			result(i, j) /= diagonal;
		}

		// Eliminate all other entries in this column
		for (int row = 0; row < n; ++row) {
			if (row != i) {
				const double factor = result(row, i);
				for (int col = 0; col < n * 2; ++col) {
					result(row, col) -= factor * result(i, col);
				}
			}
		}
	}

	// Extract the right half as the inverse
	matrix inv(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			inv(i, j) = result(i, n + j);
		}
	}

	return inv;
}

void matrix::fill_five_diagonals(const double first, const double second, const double third) {
	const int range = std::min(rows_, cols_);
	for (int i = 0; i < range; ++i) {
		// A1
		(*this)(i, i) = first;

		// A2
		if (i < range - 1)
			(*this)(i, i + 1) = second;

		if (i > 0)
			(*this)(i, i - 1) = second;

		// A3
		if (i < range - 2)
			(*this)(i, i + 2) = third;

		if (i > 1)
			(*this)(i, i - 2) = third;
	}
}


std::ostream& operator<<(std::ostream& os, const matrix& m) {
	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.cols(); ++j) {
			os << m(i, j) << ' ';
		}
		os << '\n';
	}
	return os;
}
