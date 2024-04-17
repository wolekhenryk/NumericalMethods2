#include "gauss_seidel.h"

#include <chrono>

std::tuple<matrix, int, int, std::vector<double>> gauss_seidel::solve(const matrix& a, const matrix& b) {
	const auto n = a.rows();
	const auto start = std::chrono::high_resolution_clock::now();

	matrix x(n, 1);
	x.fill(1); // Initial guess: All entries set to 1

	std::vector<double> norms;

	int iteration = 0;
	for (; iteration < MAX_ITERATIONS; iteration++) {
		matrix old_x = x; // Copy current x for convergence check

		for (int i = 0; i < n; i++) {
			double sum = b(i, 0); // Start with the corresponding b element
			for (int j = 0; j < n; j++) {
				if (i != j) {
					sum -= a(i, j) * x(j, 0); // Use the most recently updated values
				}
			}
			x(i, 0) = sum / a(i, i); // Update x(i) immediately for use in next row calculations
		}

		// Compute the norm of the difference vector between successive iterations
		matrix diff = x - old_x;
		double error = diff.norm();
		norms.push_back(error);

		if (error < EPSILON)
			break;
	}

	const auto end = std::chrono::high_resolution_clock::now();
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	return std::make_tuple(x, iteration, duration, norms);
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
