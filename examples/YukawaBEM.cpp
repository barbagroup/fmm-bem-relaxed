/** @file serialBEM.cpp
 * @brief Testing and debugging script for FMM-BEM applications
 */


#include "FMM_plan.hpp"
#include "YukawaCartesianBEM.hpp"
#include "Triangulation.hpp"
#include "gmres.hpp"

#include <sys/time.h>

double get_time()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)(tv.tv_sec + 1e-6*tv.tv_usec);
}

struct ProblemOptions
{
  typedef enum { PHI_SET, DPHIDN_SET } BoundaryCondition;
  double kappa_ = 0.;
  BoundaryCondition bc_ = PHI_SET;
  double value_ = 1.;
  int recursions = 4;

  ProblemOptions() : kappa_(0.125), bc_(PHI_SET), value_(1.) {};
  ProblemOptions(double kappa) : kappa_(kappa), bc_(PHI_SET), value_(1.) {};
  ProblemOptions(double kappa, double value) : kappa_(kappa), bc_(PHI_SET), value_(value) {};
  ProblemOptions(double kappa, BoundaryCondition bc, double value) : kappa_(kappa), bc_(bc), value_(value) {};

  double getKappa() { return kappa_; };
  double getValue() { return value_; };
  BoundaryCondition getBC() { return bc_; };
};

void printHelpAndExit()
{
  printf("serialBEM : FMM-BEM for Potential problems\n");
  printf("\nUsage: ./serialBEM <options>\n\n");
  printf("Options:\n");
  printf("-kappa <double> : Set shielding constant for the Yukawa potential\n");
  printf("-theta <double> : Set MAC theta for treecode evaluators\n");
  printf("-p <double> : Number of terms in the Multipole / Local expansions\n");
  printf("-k {1,3,4,7} : Number of Gauss integration points used per panel\n");
  printf("-lazy_eval : enable 'lazy' evaluator\n");
  printf("-ncrit <int> : Maximum # of particles per Octree box\n");
  printf("-recursions <int> : number of recursive subdivisions to create a sphere - # panels = 2*4^recursions\n");
  printf("-help : print this message\n");
  std::exit(0);
}

template <typename SourceType, typename ChargeType>
void initialiseSphere(std::vector<SourceType>& panels,
                      std::vector<ChargeType>&  charges,
                      unsigned recursions = 4)
{
  (void) charges;
  create_unit_sphere(panels, recursions);
}

int main(int argc, char **argv)
{
  int numPanels= 1000, recursions = 4, p = 5, k = 3;
  double kappa = 0.;
  FMMOptions opts;
  opts.set_mac_theta(0.5);    // Multipole acceptance criteria
  opts.set_max_per_box(10);
  SolverOptions solver_options;
  bool second_kind = false;

  // parse command line args
  // check if no arguments given
  if (argc == 1) printHelpAndExit();
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-kappa") == 0 ) {
      i++;
      kappa = atof(argv[i]);
    } else if (strcmp(argv[i],"-theta") == 0) {
      i++;
      opts.set_mac_theta((double)atof(argv[i]));
    } else if (strcmp(argv[i],"-eval") == 0) {
      i++;
      if (strcmp(argv[i],"FMM") == 0) {
        opts.evaluator = FMMOptions::FMM;
      } else if (strcmp(argv[i],"TREE") == 0) {
        opts.evaluator = FMMOptions::TREECODE;
      } else {
        printf("[W]: Unknown evaluator type: \"%s\"\n",argv[i]);
      }
    } else if (strcmp(argv[i],"-p") == 0) {
      i++;
      p = atoi(argv[i]);
    } else if (strcmp(argv[i],"-k") == 0) {
      i++;
      k = atoi(argv[i]);
    } else if (strcmp(argv[i],"-lazy_eval") == 0) {
      opts.lazy_evaluation = true;
    } else if (strcmp(argv[i],"-ncrit") == 0) {
      i++;
      opts.set_max_per_box((unsigned)atoi(argv[i]));
    } else if (strcmp(argv[i],"-recursions") == 0) {
      i++;
      recursions = atoi(argv[i]);
    } else if (strcmp(argv[i],"-second_kind") == 0) {
      second_kind = true;
    } else if (strcmp(argv[i],"-fixed_p") == 0) {
      solver_options.variable_p = false;
    } else if (strcmp(argv[i],"-help") == 0) {
      printHelpAndExit();
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
      printHelpAndExit();
    }
  }

  double tic, toc;
  tic = get_time();
  // Init the FMM Kernel
  typedef YukawaCartesianBEM kernel_type;
  kernel_type K(p, kappa, k);

  // useful typedefs
  typedef kernel_type::point_type  point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  // Init points and charges
  std::vector<source_type> panels(numPanels);
  std::vector<charge_type> charges(numPanels);
  initialiseSphere(panels, charges, recursions); //, ProblemOptions());

  // run case solving for Phi (instead of dPhi/dn) // Second-kind equation
  if (second_kind)
    for (auto& it : panels) it.switch_BC();

  // set constant Phi || dPhi/dn for each panel
  charges.resize(panels.size());
  charges = std::vector<charge_type>(panels.size(),1.);

  // Build the FMM structure
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, panels, opts);

  // generate the RHS and initial condition
  std::vector<charge_type> x(panels.size(),0.);

  // generate RHS using temporary FMM plan
  std::vector<result_type> b(panels.size(),0.);
  {
    for (auto& it : panels) it.switch_BC();
    FMM_plan<kernel_type> rhs_plan = FMM_plan<kernel_type>(K,panels,opts);
    b = rhs_plan.execute(charges,p);
    for (auto& it : panels) it.switch_BC();
  }

  toc = get_time();
  double setup_time = toc-tic;

  // Solve the system using GMRES
  // generate the Preconditioner
  tic = get_time();
  Preconditioners::Diagonal<charge_type> M(K,panels.begin(),panels.end());
  fmm_gmres(plan, x, b, solver_options, M);
  // direct_gmres(K, panels, x, b, SolverOptions());
  toc = get_time();
  double solve_time = toc-tic;

  printf("\nTIMING:\n");
  printf("\tsetup : %.4es\n",setup_time);
  printf("\tsolve : %.4es\n",solve_time);

  // check errors -- analytical solution for dPhi/dn = 1.
  double e = 0.;
  double e2 = 0.;
  double an = 1.;
  for (auto xi : x) { e += (xi-an)*(xi-an); e2 += an*an; }

  printf("error: %.3e\n",sqrt(e/e2));
}

