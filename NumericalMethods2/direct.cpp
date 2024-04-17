#include "direct.h"

#include <chrono>

std::tuple<matrix, int, int, std::vector<double>> direct::solve(const matrix& a, const matrix& b) {
	const auto start = std::chrono::high_resolution_clock::now();

	auto [l, u] = lu_decomposition(a);
	const auto y = forward_substitution(l, b);
	auto x = backward_substitution(u, y);

	const auto norm = (a * x - b).norm();

	const auto end = std::chrono::high_resolution_clock::now();
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	return std::make_tuple(x, 0, duration, std::vector{norm});
}

std::tuple<matrix, matrix> direct::lu_decomposition(const matrix& a) {
	const auto n = a.cols();

	auto l = matrix::identity(n);
	auto u = a;

	for (int i = 0; i < n; i++) {

#pragma omp parallel for
		for (int j = i + 1; j < n; j++) {
			if (u(i, i) == 0.0f)
				throw std::runtime_error("Division by zero");

			const auto factor = u(j, i) / u(i, i);
			l(j, i) = factor;

			for (int k = i; k < n; k++)
				u(j, k) -= factor * u(i, k);
		}
	}

	return std::make_tuple(l, u);
}

matrix direct::forward_substitution(const matrix& l, const matrix& b) {
	const auto n = l.cols();

	matrix y(n, 1);
	y.fill(0);

	for (int i = 0; i < n; i++) {
		y(i, 0) = b(i, 0);

		for (int j = 0; j < i; j++)
			y(i, 0) -= l(i, j) * y(j, 0);

		y(i, 0) /= l(i, i);
	}

	return y;
}

matrix direct::backward_substitution(const matrix& u, const matrix& y) {
	const auto n = u.cols();

	matrix x(n, 1);
	x.fill(0);

	for (int i = n - 1; i >= 0; i--) {
		x(i, 0) = y(i, 0);
		for (int j = i + 1; j < n; j++)
			x(i, 0) -= u(i, j) * x(j, 0);

		x(i, 0) /= u(i, i);
	}

	return x;
}
