#pragma once
#include "solution.h"

class jacobi final : public solution {
public:
	std::tuple<matrix, int, int, std::vector<double>> solve(const matrix& a, const matrix& b) override;
};
