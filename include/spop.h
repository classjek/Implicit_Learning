#ifndef SPOP_H
#define SPOP_H

#include <string>
#include <vector>
#include <tuple>

#include "metrics.h"
#include "domain.h"

// Simple wrapper to call SparsePOP
void solveWithSparsePOP(std::string& gmsFilePath, std::tuple<int,int, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<double>, std::vector<domain::BoundConstraint>>& fromGen, metrics::Checkpoint& cp);

#endif