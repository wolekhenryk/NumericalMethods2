#include <chrono>
#include <iostream>
#include <windows.h>

#include "matrix.h"

constexpr int MATRIX_SIZE = 999;
constexpr int INDEX_NUMBER = 193399;

constexpr int E_PARAM = INDEX_NUMBER / 1000 % 10;
constexpr int F_PARAM = INDEX_NUMBER / 10000 % 10;

void fill_vector_matrix(matrix&);

matrix solve_system(const matrix& a, const matrix& b);

int main() {
	SetConsoleOutputCP(CP_UTF8);

	const auto start = std::chrono::high_resolution_clock::now();

	constexpr int a1 = 5 + E_PARAM;
	constexpr int a2 = -1;
	constexpr int a3 = -1;

	matrix a(MATRIX_SIZE, MATRIX_SIZE);
	a.fill_five_diagonals(a1, a2, a3);

	matrix b(MATRIX_SIZE, 1);
	fill_vector_matrix(b);

	const auto solution = solve_system(a, b);
	std::cout << solution;

	const auto end = std::chrono::high_resolution_clock::now();

	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Time: " << duration.count() << "ms\n";

	return 0;
}

// Zadanie A
void fill_vector_matrix(matrix& m) {
	for (int i = 1; i <= m.rows(); ++i) {
		m(i - 1, 0) = std::sin(i * (F_PARAM + 1));
	}
}

matrix solve_system(const matrix& a, const matrix& b) {
	matrix result = a.inverse() * b;
	return result;
}
