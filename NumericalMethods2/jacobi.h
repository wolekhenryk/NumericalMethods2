#pragma once
#include "solution.h"

class jacobi final : public solution {
public:
	std::tuple<matrix, int, int> solve(const matrix& a, const matrix& b) override;
};
