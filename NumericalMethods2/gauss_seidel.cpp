#include "gauss_seidel.h"

#include <chrono>

std::tuple<matrix, int, int> gauss_seidel::solve(const matrix& a, const matrix& b) {
	const auto n = a.rows();

	const auto start = std::chrono::high_resolution_clock::now();

	matrix x(n, 1);
	x.fill(1);

	const auto l = get_lower_triangle(a);
	const auto u = get_upper_triangle(a);
	const auto d = get_diagonal(a);

	const auto m = (d - l).inverse();

	int iteration = 0;
	for (; iteration < MAX_ITERATIONS; iteration++) {
		for (int i = 0; i < n; i++) {
			double sum = b(i, 0);
			for (int j = 0; j < n; j++) {
				if (i != j)
					sum -= a(i, j) * x(j, 0);
			}
			x(i, 0) = sum / a(i, i);
		}

		if ((a * x - b).norm() < EPSILON)
			break;
	}

	const auto end = std::chrono::high_resolution_clock::now();
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	return {x, iteration, duration};
}

matrix gauss_seidel::get_lower_triangle(const matrix& a) {
	const auto n = a.rows();
	matrix lower_triangle(n, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			lower_triangle(i, j) = a(i, j);
		}
	}

	return lower_triangle;
}

matrix gauss_seidel::get_upper_triangle(const matrix& a) {
	const auto n = a.rows();
	matrix upper_triangle(n, n);

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			upper_triangle(i, j) = a(i, j);
		}
	}

	return upper_triangle;
}


matrix gauss_seidel::get_diagonal(const matrix& a) {
	const auto n = a.rows();
	matrix diagonal(n, n);

	for (int i = 0; i < n; i++) {
		diagonal(i, i) = a(i, i);
	}

	return diagonal;
}
