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
		double max_change = 0.0;

		for (int i = 0; i < n; i++) {
			const double old_value = x(i, 0); // Store old value of x(i) for change calculation
			double sum = b(i, 0); // Start with the corresponding b element
			for (int j = 0; j < n; j++) {
				if (i != j) {
					sum -= a(i, j) * x(j, 0); // Use the most recently updated values
				}
			}
			x(i, 0) = sum / a(i, i); // Update x(i) immediately for use in next row calculations
			if (const double change = abs(x(i, 0) - old_value); change > max_change) {
				max_change = change;
			}
		}

		norms.push_back(max_change); // Push back the max change instead of the norm

		if (max_change < EPSILON)
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
