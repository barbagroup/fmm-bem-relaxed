/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <FMM_plan.hpp>
#include <SphericalLaplaceKernel.hpp>
#include <UnitKernel.hpp>
// #include <CartesianLaplaceKernel.hpp>
#include <CartesianLaplaceKernel2.hpp>
#include <CartesianYukawaKernel.hpp>

// modify error checking for counting kernel
// TODO: Do this much better...
// #define UNIT_KERNEL
//#define SPH_KERNEL
//#define CART_KERNEL
#define YUKAWA_KERNEL

template <typename Box>
void print_box(const Box& b, std::string padding = std::string()) {
  std::cout << padding << "Box " << b.index()
            << " (Level " << b.level() << ", Parent " << b.parent().index() << "): "
            << b.morton_index() << "    " << b.center() << "\n";

  padding.append(2,' ');
  if (!b.is_leaf()) {
    for (auto ci = b.child_begin(); ci != b.child_end(); ++ci)
      print_box(*ci, padding);
  } else {
    for (auto ci = b.body_begin(); ci != b.body_end(); ++ci)
      std::cout << padding << "Point " << ci->index() << ": "
		<< ci->morton_index() << "\t" << ci->point() << "\n";
  }
}

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}

int main(int argc, char **argv)
{
  int numBodies = 100;
  bool checkErrors = true;
  bool printBox = false;
  FMMOptions opts;
  opts.set_theta(1 / sqrtf(4));    // Multipole acceptance criteria
  opts.NCRIT = 10;

  // parse command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numBodies = atoi(argv[i]);
    } else if (strcmp(argv[i],"-theta") == 0) {
      i++;
      opts.THETA = (double)atof(argv[i]);
      opts.set_theta(opts.THETA);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    } else if (strcmp(argv[i],"-bottomup") == 0) {
      opts.tree = FMMOptions::BOTTOMUP;
    } else if (strcmp(argv[i],"-evaluator") == 0) {
      i++;
      if (strcmp(argv[i],"FMM") == 0) {
        opts.evaluator = FMMOptions::FMM;
      } else if (strcmp(argv[i],"TREECODE") == 0) {
        opts.evaluator = FMMOptions::TREECODE;
      } else {
        printf("[W]: Unknown evaluator type: \"%s\"\n",argv[i]);
      }
    } else if (strcmp(argv[i],"-ncrit") == 0) {
      i++;
      opts.NCRIT = (unsigned)atoi(argv[i]);
    } else if (strcmp(argv[i],"-printbox") == 0) {
      printBox = true;
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
    }
  }

  // Init the FMM Kernel
#ifdef SPH_KERNEL
  typedef SphericalLaplaceKernel kernel_type;
  kernel_type K(5);
#endif
#ifdef CART_KERNEL
  typedef CartesianLaplaceKernel<5> kernel_type;
  kernel_type K;
#endif
#ifdef YUKAWA_KERNEL
  typedef CartesianYukawaKernel kernel_type;
  kernel_type K(6,0);
#endif
#ifdef UNIT_KERNEL
  typedef UnitKernel kernel_type;
  kernel_type K;
#endif
  typedef kernel_type::point_type point_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;


  // Init points and charges
  std::vector<point_type> points(numBodies);
  for (int k = 0; k < numBodies; ++k)
    points[k] = point_type(drand(), drand(), drand());
  //std::vector<point_type> target_points = points;

  std::vector<charge_type> charges(numBodies);
  for (int k = 0; k < numBodies; ++k)
    charges[k] = drand();

  // Build the FMM
  //fmm_plan plan = fmm_plan(K, bodies, opts);
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, points, opts);
  if (printBox) {
    print_box(plan.otree.root());
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
