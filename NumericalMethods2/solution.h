#pragma once
#include "matrix.h"

constexpr double EPSILON = 1e-9;
constexpr int MAX_ITERATIONS = 1e3;

class solution {
public:
	virtual std::tuple<matrix, int, int, std::vector<double>> solve(const matrix& a, const matrix& b) = 0;
	virtual ~solution() = default;
};
