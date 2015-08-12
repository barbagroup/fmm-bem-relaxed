/** @file serialun.cpp
 * @brief Testing and debugging script
 */

//#define DEBUG

#include <FMM_plan.hpp>
//#include <LaplaceSphericalModified.hpp>
#include <LaplaceSpherical.hpp>
//#include <UnitKernel.hpp>
//#include <KernelSkeleton.hpp>
//#include <KernelSkeletonMixed.hpp>
#include <LaplaceCartesian.hpp>
//#include <YukawaCartesian.hpp>
#include <YukawaSpherical.hpp>
#include "StokesSpherical.hpp"

#include <cstring>

// modify error checking for counting kernel
// TODO: Do this much better...
//#define SKELETON_KERNEL
//#define UNIT_KERNEL
// #define SPH_KERNEL
//#define CART_KERNEL
//#define YUKAWA_KERNEL
#define YUKAWA_SPH
//#define STOKES_SPH

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}

// Random number in [A,B)
inline double drand(double A, double B) {
  return A + (B-A) * drand();
}

int main(int argc, char **argv)
{
  //std::vector<std::string> args(argv+1, argv+argc);
  //std::cout << args[0] << '\n' << args[1] << '\n' << args[2] << std::endl;

  int numBodies = 1000, p=5;
  bool checkErrors = true;
  double beta = 0.125;

  FMMOptions opts = get_options(argc, argv);

  //opts.local_evaluation = true;
  //opts.sparse_local = true;

  // parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numBodies = atoi(argv[i]);
    } else if (strcmp(argv[i],"-p") == 0) {
      i++;
      p = atoi(argv[i]);
    } else if (strcmp(argv[i],"-beta") == 0) {
      i++;
      beta = atof(argv[i]);
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
  kernel_type K(p);
#endif
#ifdef CART_KERNEL
  typedef LaplaceCartesian<5> kernel_type;
  kernel_type K;
#endif
#ifdef YUKAWA_KERNEL
  typedef YukawaCartesian kernel_type;
  kernel_type K(p,beta);
#endif
#ifdef YUKAWA_SPH
  typedef YukawaSpherical kernel_type;
  kernel_type K(p,beta);
#endif
#ifdef UNIT_KERNEL
  typedef UnitKernel kernel_type;
  kernel_type K;
#endif
#ifdef STOKES_SPH
  typedef StokesSpherical kernel_type;
  kernel_type K(p);
#endif
  // if not using a Yukawa kernel, quiet warnings
#if !defined(YUKAWA_KERNEL) && !defined(YUKAWA_SPH)
  (void) beta;
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
  for (int k = 0; k < numBodies; ++k) {
#if defined(STOKES_SPH)
    charges[k] = charge_type(1,1,1); // charge_type(drand(),drand(),drand());
#else
    charges[k] = drand();
#endif
  }

  // Build the FMM
  //fmm_plan plan = fmm_plan(K, bodies, opts);
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, points, opts);

  // Execute the FMM
  //fmm_execute(plan, charges, target_points);
  std::vector<result_type> result = plan.execute(charges);

  // Check the result
  // TODO: More elegant
  if (checkErrors) {
    // choose a number of samples to use
    int numSamples = std::min(numBodies, 1000);
    std::vector<int> sample_map(numSamples);
    std::vector<point_type> sample_targets(numSamples);
    std::vector<result_type> exact(numSamples);

    // sample the space (needs to be done better)
    for (int i=0; i<numSamples; i++) {
      int sample = i >> 15 % numBodies;
      sample_map[i] = sample;
      sample_targets[i] = points[sample];
    }

    // Compute the result with a direct matrix-vector multiplication
    Direct::matvec(K, points.begin(), points.end(), charges.begin(),
                   sample_targets.begin(), sample_targets.end(), exact.begin());

#if defined(SPH_KERNEL) || defined(CART_KERNEL) || defined(YUKAWA_KERNEL)
    result_type rdiff, rnorm;
    for (unsigned k = 0; k < exact.size(); ++k) {
      auto exact_val = exact[k];
      auto fmm_val   = result[sample_map[k]];
      // printf("[%03d] exact: %lg, FMM: %lg\n", k, exact_val[0], fmm_val[0]);

      rdiff += (fmm_val - exact_val) * (fmm_val - exact_val);
      rnorm += exact_val * exact_val;
    }

    printf("Error (pot) : %.4e\n", sqrt(rdiff[0] / rnorm[0]));
    printf("Error (acc) : %.4e\n", sqrt((rdiff[1]+rdiff[2]+rdiff[3]) /
					(rnorm[1]+rnorm[2]+rnorm[3])));
#endif
#if defined(YUKAWA_SPH)
    result_type rdiff = 0., rnorm = 0.;
    for (unsigned k = 0; k < exact.size(); ++k) {
      // printf("[%03d] exact: %lg, FMM: %lg\n", k, exact[k], result[k]);
      auto exact_val = exact[k];
      auto fmm_val   = result[sample_map[k]];

      rdiff = (fmm_val - exact_val) * (fmm_val - exact_val);
      rnorm = exact_val * exact_val;
    }

    printf("Error (pot) : %.4e\n", sqrt(rdiff / rnorm));
#endif
#if defined(STOKES_SPH)
    result_type rdiff = result_type(0.), rnorm = result_type(0.);
    for (unsigned k = 0; k < exact.size(); ++k) {
      auto exact_val = exact[k];
      auto fmm_val = result[sample_map[k]];

      for (unsigned i=0; i < 3; i++) {
        rdiff[i] += (fmm_val[i]-exact_val[i])*(fmm_val[i]-exact_val[i]);
        rnorm[i] += exact_val[i]*exact_val[i];
      }
    }
    auto div = rdiff/rnorm;
    printf("Error (u) : %.4e, (v) : %.4e, (w) : %.4e\n",sqrt(div[0]),sqrt(div[1]),sqrt(div[2]));
#endif
#ifdef UNIT_KERNEL
    int wrong_results = 0;
    for (unsigned k = 0; k < exact.size(); ++k) {
      auto exact_val = exact[k];
      auto fmm_val   = result[sample_map[k]];

      // printf("[%03d] exact: %lg, FMM: %lg\n", k, exact[k], result[k]);

      if ((exact_val - fmm_val) / exact_val > 1e-13)
	++wrong_results;
    }
    printf("Wrong counts: %d\n", wrong_results);
#endif
  }
}
