#ifndef SPOP_H
#define SPOP_H

#include <string>
#include <vector>
#include <tuple>

// Simple wrapper to call SparsePOP
void solveWithSparsePOP(std::string& gmsFilePath, std::tuple<int,int, std::vector<int>, std::vector<int>, std::vector<int>>& fromGen);

#endif