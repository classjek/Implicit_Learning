#include "executor.h"
#include <iostream>
#include <omp.h>

int Executor::run() {
  std::cout << "[executor] degree_d=" << cfg_.degree_d
            << " eps_bisect=" << cfg_.eps_bisect
            << " omp_threads=" << cfg_.omp_threads
            << " solver_threads=" << cfg_.solver_threads
            << std::endl;

  omp_set_num_threads(cfg_.omp_threads);

  std::cout << "[executor] omp_get_num_procs() = " << omp_get_num_procs() << '\n'
            << "[executor] omp_get_max_threads() = " << omp_get_max_threads() << '\n';

  // One thread here for each observation
  // for each observation, generate bounds for each equivalence class, assign to corresponding constraint
  #pragma omp parallel
  {
    #pragma omp critical 
    {
      std::cout << omp_get_thread_num() << " of " << omp_get_num_threads() << " threads\n";
    }
    // #pragma omp single
    // std::cout << "[executor] OpenMP using " << omp_get_num_threads() << " threads\n";
  }

  int m = 8;  // placeholder
  #pragma omp parallel for schedule(dynamic,1)
  for (int j = 0; j < m; ++j) {
    // per-sample work
  }

  return 0;
}