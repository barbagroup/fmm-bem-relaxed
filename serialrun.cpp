/** @file serialun.cpp
 * @brief Testing and debugging script
 */

//#define DEBUG

#include <FMM_plan.hpp>
//#include <LaplaceSphericalModified.hpp>
#include <LaplaceSpherical.hpp>
#include <UnitKernel.hpp>
#include <KernelSkeleton.hpp>
//#include <KernelSkeletonMixed.hpp>
#include <LaplaceCartesian.hpp>
#include <YukawaCartesian.hpp>

#include <string.h>

// modify error checking for counting kernel
// TODO: Do this much better...
//#define SKELETON_KERNEL
//#define UNIT_KERNEL
#define SPH_KERNEL
//#define CART_KERNEL
//#define YUKAWA_KERNEL

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

  FMMOptions opts = get_options(argc, argv);

  // parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numBodies = atoi(argv[i]);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    }
  }

  // Init the FMM Kernel
#ifdef SKELETON_KERNEL
  typedef KernelSkeleton kernel_type;
  kernel_type K;
#endif
#ifdef SPH_KERNEL
  typedef LaplaceSpherical kernel_type;
  kernel_type K(5);
#endif
#ifdef CART_KERNEL
  typedef LaplaceCartesian<5> kernel_type;
  kernel_type K;
#endif
#ifdef YUKAWA_KERNEL
  typedef YukawaCartesian kernel_type;
  kernel_type K(6,0.5);
#endif
#ifdef UNIT_KERNEL
  typedef UnitKernel kernel_type;
  kernel_type K;
#endif
  typedef kernel_type::point_type point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;


  // Init points and charges
  std::vector<source_type> points(numBodies);
  for (int k = 0; k < numBodies; ++k)
    points[k] = source_type(drand(), drand(), drand());
  //std::vector<point_type> target_points = points;

  std::vector<charge_type> charges(numBodies);
  for (int k = 0; k < numBodies; ++k)
    charges[k] = drand();

  // Build the FMM
  //fmm_plan plan = fmm_plan(K, bodies, opts);
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, points, opts);

  // Execute the FMM
  //fmm_execute(plan, charges, target_points);
  std::vector<result_type> result = plan.execute(charges);

  // Check the result
  // TODO: More elegant
  if (checkErrors) {
    std::vector<result_type> exact(numBodies);

    // Compute the result with a direct matrix-vector multiplication
    Direct::matvec(K, points, charges, exact);

#if defined(SPH_KERNEL) || defined(CART_KERNEL) || defined(YUKAWA_KERNEL)
    result_type rdiff, rnorm;
    for (unsigned k = 0; k < result.size(); ++k) {
      printf("[%03d] exact: %lg, FMM: %lg\n", k, exact[k][0], result[k][0]);

      rdiff = (result[k] - exact[k]) * (result[k] - exact[k]);
      rnorm = exact[k] * exact[k];
    }

    printf("Error (pot) : %.4e\n", sqrt(rdiff[0] / rnorm[0]));
    printf("Error (acc) : %.4e\n", sqrt((rdiff[1]+rdiff[2]+rdiff[3]) /
					(rnorm[1]+rnorm[2]+rnorm[3])));
#endif
#ifdef UNIT_KERNEL
    int wrong_results = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
      printf("[%03d] exact: %lg, FMM: %lg\n", k, exact[k], result[k]);

      if ((exact[k] - result[k]) / exact[k] > 1e-13)
	++wrong_results;
    }
    printf("Wrong counts: %d\n", wrong_results);
#endif
  }
}
