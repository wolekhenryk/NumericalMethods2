#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <windows.h>

#include "direct.h"
#include "gauss_seidel.h"
#include "jacobi.h"
#include "matrix.h"

constexpr int INDEX_NUMBER = 193399;
constexpr int MATRIX_SIZE = 900 + INDEX_NUMBER % 100;

constexpr int E = INDEX_NUMBER / 100 % 10;
constexpr int F = INDEX_NUMBER / 1000 % 10;

constexpr int DISPLAY_PRECISION = 5;

constexpr std::array MATRIX_SIZES = { 100, 500, 1000, 2000, 3000, 6000 };

void present_jacobi_method(const matrix& a, const matrix& b);

void present_gauss_method(const matrix& a, const matrix& b);

void present_lu_method(const matrix& a, const matrix& b);

int main() {
	SetConsoleOutputCP(CP_UTF8);

	constexpr int a1 = 5 + E;
	constexpr int a2 = -1;
	constexpr int a3 = -1;

	matrix a(MATRIX_SIZE, MATRIX_SIZE);
	a.fill_five_diagonals(a1, a2, a3);

	matrix b(MATRIX_SIZE, 1);
	b.fill_vector_with_sine(F);

	//for (const auto size : MATRIX_SIZES) {
	//	matrix a(size, size);
	//	a.fill_five_diagonals(a1, a2, a3);

	//	matrix b(size, 1);
	//	b.fill_vector_with_sine(F_PARAM);

	//	std::cout << size << '\t';

		//present_jacobi_method(a, b);
		//present_gauss_method(a, b);
	//	present_lu_method(a, b);

	//	std::cout << '\n';
	//}

	present_jacobi_method(a, b);
	present_gauss_method(a, b);
	//present_lu_method(a, b);

	//std::cout << "N" << '\t' << "Jacobi time, norm" << '\t';
	//std::cout << "Gauss time, norm" << '\t';
	//std::cout << "LU time, norm" << '\n';

	//for (const auto size : MATRIX_SIZES) {
	//	matrix a(size, size);
	//	a.fill_five_diagonals(a1, a2, a3);

	//	matrix b(size, 1);
	//	b.fill_vector_with_sine(F);

	//	std::cout << size << '\t';
	//	present_jacobi_method(a, b);
	//	present_gauss_method(a, b);
	//	present_lu_method(a, b);
	//	std::cout << '\n';
	//}

	return 0;
}

void present_jacobi_method(const matrix& a, const matrix& b) {
	//std::cout << "\n\nJacobi Method\n\n";

	const std::unique_ptr<solution> jacobi_solver = std::make_unique<jacobi>();
	const auto [sol_jacobi, iterations, duration, norms] = jacobi_solver->solve(a, b);

	//for (int i = 0; i < sol_jacobi.rows(); i++) {
	//			std::cout << sol_jacobi(i, 0) << '\n';
	//}

	//for (const auto norm : norms)
	//	std::cout << norm << '\n';

	std::cout << "\n\n";

	std::cout << "Iterations: " << iterations
		<< "\nDuration: " << duration
		<< "ms\nLast norm: " << std::scientific << std::setprecision(DISPLAY_PRECISION) << norms.back() << '\n';
}

void present_gauss_method(const matrix& a, const matrix& b) {
	//std::cout << "\n\nGauss method\n\n";

	const std::unique_ptr<solution> gauss_solver = std::make_unique<gauss_seidel>();
	const auto [sol_gauss, iterations, duration, norms] = gauss_solver->solve(a, b);

	//for (const auto norm : norms)
	//	std::cout << norm << '\n';

	std::cout << "\n\n";

	std::cout << "Iterations: " << iterations
		<< "\nDuration: " << duration
		<< "ms\nLast norm: " << std::scientific << std::setprecision(DISPLAY_PRECISION) << norms.back() << '\n';

	//std::cout << duration << '\t' << iterations << '\n';
}

void present_lu_method(const matrix& a, const matrix& b) {
	//std::cout << "\n\nLU method\n\n";

	const std::unique_ptr<solution> lu_solver = std::make_unique<direct>();
	const auto [sol_lu, iterations, duration, norms] = lu_solver->solve(a, b);

	//std::cout << "\nDuration: " << duration
	//	<< "ms\nNorm: " << std::scientific << std::setprecision(DISPLAY_PRECISION) << norms.back() << '\n';

	std::cout << duration << '\t';
}