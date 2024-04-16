#include "jacobi.h"

#include <chrono>

std::tuple<matrix, int, int> jacobi::solve(const matrix& a, const matrix& b) {
	const auto start = std::chrono::high_resolution_clock::now();

	const int n = a.rows();
	matrix x(n, 1);
	matrix x_new(n, 1);

	int iteration = 0;
	for (; iteration < MAX_ITERATIONS; ++iteration) {
		for (int i = 0; i < n; ++i) {
			double sum = 0.0;
			for (int j = 0; j < n; ++j) {
				if (i != j) {
					sum += a(i, j) * x(j, 0);
				}
			}
			x_new(i, 0) = (b(i, 0) - sum) / a(i, i);
		}

		if (const double error = (x_new - x).norm(); error < EPSILON)
			break;

		x = x_new;
	}

	const auto end = std::chrono::high_resolution_clock::now();
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	return {x_new, iteration, duration};
}
