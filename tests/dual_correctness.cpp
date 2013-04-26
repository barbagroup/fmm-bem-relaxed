/** @file correctness.cpp
 * @brief Test the tree and tree traversal by running an instance
 * of the UnitKernel with random points and charges
 */

#include "FMM_plan.hpp"
#include "UnitKernel.hpp"

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}

// Random number in [A,B)
inline double drand(double A, double B) {
  return (B-A) * drand() + A;
}

int main(int argc, char **argv)
{
  int numBodies = 1000;
  bool checkErrors = true;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numBodies = atoi(argv[i]);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    }
  }

  // Init the FMM Kernel and options
  FMMOptions opts = get_options(argc, argv);
  typedef UnitKernel kernel_type;
  kernel_type K;

  typedef kernel_type::point_type point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  // Init sources
  std::vector<source_type> sources(numBodies);
  for (int k = 0; k < numBodies; ++k)
    sources[k] = source_type(drand(), drand(), drand());
  // Init targets
  std::vector<target_type> targets(numBodies);
  for (int k = 0; k < numBodies; ++k)
    targets[k] = target_type(drand(), drand(), drand());
  // Init charges
  std::vector<charge_type> charges(numBodies);
  for (int k = 0; k < numBodies; ++k)
    charges[k] = drand();

  // Build the FMM
  //fmm_plan plan = fmm_plan(K, bodies, opts);
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, sources, targets, opts);

  // Execute the FMM
  //fmm_execute(plan, charges, target_points);
  std::vector<result_type> result = plan.execute(charges);

  // Check the result
  if (checkErrors) {
    std::vector<result_type> exact(numBodies);

    // Compute the result with a direct matrix-vector multiplication
    Direct::matvec(K, sources, charges, targets, exact);

    int wrong_results = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
      printf("[%03d] exact: %lg, FMM: %lg\n", k, exact[k], result[k]);

      if ((exact[k] - result[k]) / exact[k] > 1e-13)
        ++wrong_results;
    }
    printf("Wrong counts: %d\n", wrong_results);
  }
}
