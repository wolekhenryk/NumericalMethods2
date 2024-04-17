#pragma once

#include "solution.h"

class direct final : public solution {
public:
	std::tuple<matrix, int, int, std::vector<double>> solve(const matrix& a, const matrix& b) override;

private:
	static std::tuple<matrix, matrix> lu_decomposition(const matrix& a);

	static matrix forward_substitution(const matrix& l, const matrix& b);

	static matrix backward_substitution(const matrix& u, const matrix& y);
};
