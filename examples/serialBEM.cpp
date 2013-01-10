/** @file serialBEM.cpp
 * @brief Testing and debugging script for FMM-BEM applications
 */


#include "FMM_plan.hpp"
#include "LaplaceSphericalBEM.hpp"
#include "Triangulation.hpp"
#include "gmres.hpp"
// #include "Preconditioner.hpp"

struct SolverOptions
{
  double residual;
  int max_iters, restart;

  SolverOptions() : residual(1e-5), max_iters(50), restart(50) {};
};

struct ProblemOptions
{
  typedef enum { PHI_SET, DPHIDN_SET } BoundaryCondition;
  BoundaryCondition bc_ = PHI_SET;
  double value_ = 1.;
  int recursions = 4;

  ProblemOptions() : bc_(PHI_SET), value_(1.) {};
  ProblemOptions(double value) : bc_(PHI_SET), value_(value) {};
  ProblemOptions(BoundaryCondition bc, double value) : bc_(bc), value_(value) {};

  double getValue() { return value_; };
  BoundaryCondition getBC() { return bc_; };
};

void printHelpAndExit()
{
  printf("serialBEM : FMM-BEM for Potential problems\n");
  printf("\nUsage: ./serialBEM <options>\n\n");
  printf("Options:\n");
  printf("-theta <double> : Set MAC theta for treecode evaluators\n");
  printf("-eval {FMM,TREE} : Choose either FMM or treecode evaluator\n");
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
  FMMOptions opts;
  opts.set_mac_theta(0.5);    // Multipole acceptance criteria
  opts.set_max_per_box(10);

  // parse command line args
  // check if no arguments given
  if (argc == 1) printHelpAndExit();
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numPanels = atoi(argv[i]);
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
    } else if (strcmp(argv[i],"-help") == 0) {
      printHelpAndExit();
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
      printHelpAndExit();
    }
  }

  // Init the FMM Kernel
  typedef LaplaceSphericalBEM kernel_type;
  kernel_type K(p,k);

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

  // run case solving for Phi (instead of dPhi/dn)
  for (auto& it : panels) it.switch_BC();

  // set constant Phi || dPhi/dn for each panel
  charges.resize(panels.size());
  charges = std::vector<charge_type>(panels.size(),1.);

  // Build the FMM structure
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, panels, opts);

  // generate the RHS and initial condition
  std::vector<charge_type> x(panels.size(),0.);

  // generate RHS using direct calculation
  for (auto& it : panels) it.switch_BC();
  std::vector<result_type> b(panels.size(),0.);
  Direct::matvec(K,panels,charges,b);
  for (auto& it : panels) it.switch_BC();

  // Solve the system using GMRES
  // generate the Preconditioner
  Preconditioners::Diagonal<charge_type> M(K,panels.begin(),panels.end());
  fmm_gmres(plan, x, b, SolverOptions(),M);
  //direct_gmres(K, panels, x, b, SolverOptions());

  // check errors -- analytical solution for dPhi/dn = 1.
  double e = 0.;
  double e2 = 0.;
  double an = 1.;
  for (auto xi : x) { e += (xi-an)*(xi-an); e2 += an*an; }
  
  printf("error: %.3e\n",sqrt(e/e2));
}

