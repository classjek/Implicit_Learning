#pragma once

struct Config {
  int    degree_d     = 2;
  double eps_bisect   = 1e-4;
  int    omp_threads  = 1;   // you can ignore if not using OpenMP yet
  int    solver_threads = 1; // ditto
  bool   verbose      = true;
};