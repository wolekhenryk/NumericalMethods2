#pragma once
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

class matrix {
	int rows_, cols_;
	std::vector<std::vector<double>> data_;

public:
	matrix(int rows, int cols);
	explicit matrix(const std::vector<std::vector<double>>& data);

	// Copy constructor
	matrix(const matrix& other);

	// Move constructor
	matrix(matrix&& other) noexcept;

	// Copy assignment
	matrix& operator=(const matrix& other);

	// Getters

	[[nodiscard]] int rows() const;

	[[nodiscard]] int cols() const;

	// Access element
	double& operator()(int row, int col);

	const double& operator()(int row, int col) const;

	[[nodiscard]] matrix operator*(const matrix& other) const;

	[[nodiscard]] matrix operator+(const matrix& other) const;

	[[nodiscard]] matrix operator-(const matrix& other) const;

	[[nodiscard]] matrix operator-() const;

	[[nodiscard]] matrix transpose() const;

	[[nodiscard]] matrix inverse() const;

	[[nodiscard]] static matrix identity(int size);

	template <typename T>
	void make_random(int start, int end);

	void fill(double value);

	[[nodiscard]] double norm() const;

	// Helper functions

	void fill_five_diagonals(double first, double second, double third);

	void fill_vector_with_sine(int f);

	// Print matrix
	friend std::ostream& operator<<(std::ostream& os, const matrix& m);
};

template <typename T>
void matrix::make_random(const int start, const int end) {
	std::mt19937 rng(std::random_device{}());

	if constexpr (std::is_same_v<T, int>) {
		std::uniform_int_distribution dist(start, end);
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				data_[i][j] = dist(rng);
			}
		}
	}
	else if constexpr (std::is_same_v<T, double>) {
		std::uniform_real_distribution<> dist(start, end);
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				data_[i][j] = dist(rng);
			}
		}
	}
	else {
		static_assert(std::is_same_v<T, int> || std::is_same_v<T, double>, "Only int and double types are allowed.");
	}
}
