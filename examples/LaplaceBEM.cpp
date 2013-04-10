/** @file serialBEM.cpp
 * @brief Testing and debugging script for FMM-BEM applications
 */


#include "FMM_plan.hpp"
#include "DirectMatvec.hpp"
#include "LaplaceSphericalBEM.hpp"
#include "Triangulation.hpp"

#include "GMRES.hpp"
#include "fmgmres.hpp"
#include "LocalPC.hpp"

#include "timing.hpp"

enum SolveType { INITIAL_SOLVE = 1,
                 GMRES_DIAGONAL = 2,
                 INNER_OUTER = 4 };

enum SOLVERS { SOLVE_GMRES, SOLVE_FGMRES } ;
enum PRECONDITIONERS { IDENTITY, DIAGONAL, LOCAL };

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
  FMMOptions opts = get_options(argc,argv);
  // opts.set_max_per_box(10);
  SolverOptions solver_options;
  bool second_kind = false;

  // solve / PC settings
  SOLVERS solver = SOLVE_GMRES;
  PRECONDITIONERS pc = IDENTITY;

  // use lazy evaluator by default
  // opts.lazy_evaluation = true;

  // parse command line args
  // check if no arguments given
  if (argc == 1) printHelpAndExit();
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numPanels = atoi(argv[i]);
    } else if (strcmp(argv[i],"-p") == 0) {
      i++;
      p = atoi(argv[i]);
      solver_options.max_p = p;
    } else if (strcmp(argv[i],"-k") == 0) {
      i++;
      k = atoi(argv[i]);
    } else if (strcmp(argv[i],"-recursions") == 0) {
      i++;
      recursions = atoi(argv[i]);
    } else if (strcmp(argv[i],"-second_kind") == 0) {
      second_kind = true;
    } else if (strcmp(argv[i],"-fixed_p") == 0) {
      solver_options.variable_p = false;
    } else if (strcmp(argv[i],"-solver_tol") == 0) {
      i++;
      solver_options.residual = (double)atof(argv[i]);
    } else if (strcmp(argv[i],"-gmres") == 0) {
      solver = SOLVE_GMRES;
    } else if (strcmp(argv[i],"-local") == 0) {
      solver = SOLVE_FGMRES;
      pc = LOCAL;
    } else if (strcmp(argv[i],"-diagonal") == 0) {
      pc = DIAGONAL;
    } else if (strcmp(argv[i],"-help") == 0) {
      printHelpAndExit();
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
      // printHelpAndExit();
    }
  }

  // opts.sparse_local = true;
  double tic, toc;
  tic = get_time();
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
  if (second_kind)
    for (auto& it : panels) it.switch_BC();

  // set constant Phi || dPhi/dn for each panel
  charges.resize(panels.size());
  charges = std::vector<charge_type>(panels.size(),1.);

  // Build the FMM structure
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, panels, opts);

  // generate the RHS and initial condition
  std::vector<charge_type> x(panels.size(),0.);

  std::vector<result_type> b(panels.size(),0.);
  // generate RHS using temporary FMM plan
  {
    for (auto& it : panels) it.switch_BC();
    FMM_plan<kernel_type> rhs_plan = FMM_plan<kernel_type>(K,panels,opts);
    b = rhs_plan.execute(charges);
    for (auto& it : panels) it.switch_BC();
  }

  toc = get_time();
  double setup_time = toc-tic;

  // Solve the system using GMRES
  // generate the Preconditioner
  tic = get_time();
  Preconditioners::Diagonal<charge_type> M(K,panels.begin(),panels.end());
  //M.print();
  SolverOptions inner_options(1e-2,1,2);
  inner_options.variable_p = true;
  // Preconditioners::FMGMRES<FMM_plan<kernel_type>,Preconditioners::Diagonal<charge_type>> inner(plan, b, inner_options, M);

  // Local preconditioner
  Preconditioners::LocalInnerSolver<FMM_plan<kernel_type>, Preconditioners::Diagonal<charge_type>> local(K, panels, b);

  // Initial low accuracy solve
  //
  /*
  double tic2, toc2;
  tic2 = get_time();
  {
    printf("Initial solve starting..\n");
    // initial solve to 1e-2 accuracy, 5 iterations, P = 2
    SolverOptions initial_solve(5e-3,50,3);
    // fmm_gmres(plan, x, b, solver_options, M);
    GMRES(plan, x, b, initial_solve, M);
    printf("Initial solve finished..\n");
  }
  toc2 = get_time();
  printf("Initial solve took: %.4es\n",toc2-tic2);
  */

  // Outer GMRES solve with diagonal preconditioner & relaxation
  solver_options.residual = 1e-5;
  FGMRESContext<result_type> context(x.size(), solver_options.restart);

  if (solver == SOLVE_GMRES && pc == IDENTITY){
    printf("Solver: GMRES\nPreconditioner: Identity\n");
    // straight GMRES, no preconditioner
    // DirectMV<kernel_type> MV(K, panels, panels);
    GMRES(plan,x,b,solver_options);
  }
  else if (solver == SOLVE_GMRES && pc == DIAGONAL) {
    printf("Solver: GMRES\nPreconditioner: Diagonal\n");
    // GMRES, diagonal preconditioner
    GMRES(plan,x,b,solver_options, M, context);
  }
  else if (solver == SOLVE_FGMRES && pc == LOCAL) {
    printf("Solver: FGMRES\nPreconditioner: Local solve\n");
    // FGMRES, Local inner solver
    FGMRES(plan,x,b,solver_options, local, context);
  }
  else {
    printf("[E] no valid solver / preconditioner option chosen\n");
    exit(0);
  }

  // GMRES(MV,x,b,solver_options, M, context);
  // FGMRES(plan,x,b,solver_options, inner, context); // , context);
  // Outer/Inner FGMRES / GMRES (Diagonal)
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

