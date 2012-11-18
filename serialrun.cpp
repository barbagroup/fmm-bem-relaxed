/** @file serialun.cpp
 * @brief Testing and debugging script
 */

#include <FMM_plan.hpp>
//#include <SphericalLaplaceKernelModified.hpp>
#include <SphericalLaplaceKernel.hpp>
#include <UnitKernel.hpp>
#include <KernelSkeleton.hpp>
//#include <KernelSkeletonMixed.hpp>
#include <CartesianLaplaceKernel.hpp>
#include <CartesianYukawaKernel.hpp>

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
  bool printBox = true;
  FMMOptions opts;
  opts.set_mac_theta(0.5);    // Multipole acceptance criteria
  opts.set_max_per_box(10);

  // parse command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numBodies = atoi(argv[i]);
    } else if (strcmp(argv[i],"-theta") == 0) {
      i++;
      opts.set_mac_theta((double)atof(argv[i]));
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    } else if (strcmp(argv[i],"-eval") == 0) {
      i++;
      if (strcmp(argv[i],"FMM") == 0) {
        opts.evaluator = FMMOptions::FMM;
      } else if (strcmp(argv[i],"TREE") == 0) {
        opts.evaluator = FMMOptions::TREECODE;
      } else {
        printf("[W]: Unknown evaluator type: \"%s\"\n",argv[i]);
      }
    } else if (strcmp(argv[i],"-lazy_eval") == 0) {
      opts.lazy_evaluation = true;
    } else if (strcmp(argv[i],"-ncrit") == 0) {
      i++;
      opts.set_max_per_box((unsigned)atoi(argv[i]));
    } else if (strcmp(argv[i],"-printbox") == 0) {
      printBox = true;
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
    }
  }

  // Init the FMM Kernel
#ifdef SKELETON_KERNEL
  typedef KernelSkeleton kernel_type;
  kernel_type K;
#endif
#ifdef SPH_KERNEL
  typedef SphericalLaplaceKernel kernel_type;
  kernel_type K(5);
#endif
#ifdef CART_KERNEL
  typedef CartesianLaplaceKernel<8> kernel_type;
  kernel_type K;
#endif
#ifdef YUKAWA_KERNEL
  typedef CartesianYukawaKernel kernel_type;
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
  if (printBox) {
    //print_box(plan.otree.root());
    //std::cout << plan.otree << "\n";
  }

  // Execute the FMM
  //fmm_execute(plan, charges, target_points);
  std::vector<result_type> result = plan.execute(charges, points);

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
