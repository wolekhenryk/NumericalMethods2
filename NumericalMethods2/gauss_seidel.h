#pragma once
#include "solution.h"

class gauss_seidel final : public solution {
public:
	std::tuple<matrix, int, int, std::vector<double>> solve(const matrix& a, const matrix& b) override;

private:
	static matrix get_lower_triangle(const matrix& a);
	static matrix get_upper_triangle(const matrix& a);
	static matrix get_diagonal(const matrix& a);
};
