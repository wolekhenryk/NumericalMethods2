#include <chrono>
#include <iostream>
#include <tuple>
#include <windows.h>

#include "gauss_seidel.h"
#include "jacobi.h"
#include "matrix.h"

constexpr int MATRIX_SIZE = 999;
constexpr int INDEX_NUMBER = 193399;

constexpr int E_PARAM = INDEX_NUMBER / 1000 % 10;
constexpr int F_PARAM = INDEX_NUMBER / 10000 % 10;

void present_jacobi_method(const matrix& a, const matrix& b);

void present_gauss_method(const matrix& a, const matrix& b);

int main() {
	SetConsoleOutputCP(CP_UTF8);

	constexpr int a1 = 5 + E_PARAM;
	constexpr int a2 = -1;
	constexpr int a3 = -1;

	matrix a(MATRIX_SIZE, MATRIX_SIZE);
	a.fill_five_diagonals(a1, a2, a3);

	matrix b(MATRIX_SIZE, 1);
	b.fill_vector_with_sine(F_PARAM);

	present_jacobi_method(a, b);
	present_gauss_method(a, b);

	return 0;
}

void present_jacobi_method(const matrix& a, const matrix& b) {
	std::cout << "\n\nJacobi Method\n\n";

	const std::unique_ptr<solution> jacobi_solver = std::make_unique<jacobi>();
	const auto [sol_jacobi, iterations, duration] = jacobi_solver->solve(a, b);

	std::cout << "Iterations: " << iterations << "\nDuration: " << duration << "ms\n";
}

void present_gauss_method(const matrix& a, const matrix& b) {
	std::cout << "\n\nGauss method\n\n";

	const std::unique_ptr<solution> gauss_solver = std::make_unique<gauss_seidel>();
	const auto [sol_gauss, iterations, duration] = gauss_solver->solve(a, b);

	std::cout << "Iterations: " << iterations << "\nDuration: " << duration << "ms\n";
}
