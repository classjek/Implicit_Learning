#pragma once
#include "config.h"

class Executor {
public:
  explicit Executor(const Config& cfg) : cfg_(cfg) {}
  // Return 0 = success for now
  int run();
private:
  Config cfg_;
};