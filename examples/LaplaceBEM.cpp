/** @file serialBEM.cpp
 * @brief Testing and debugging script for FMM-BEM applications
 */

#include "StokesSphericalBEM.hpp"

#include "FMM_plan.hpp"
#include "DirectMatvec.hpp"
#include "LaplaceSphericalBEM.hpp"
#include "Triangulation.hpp"

#include "GMRES.hpp"
#include "fmgmres.hpp"
#include "LocalPC.hpp"
#include "BlockDiagonalPC.hpp"

#include "MshReader.hpp"

#include "timing.hpp"

// simple or more complicated test
// #define BEMCPP_TEST

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
  int numPanels= 1000, recursions = 4, p = 5, k = 3, max_iterations = 500;
  FMMOptions opts = get_options(argc,argv);
  opts.sparse_local = true;
  SolverOptions solver_options;
  bool second_kind = false;
  char *mesh_name;
  bool mesh = false;

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
    } else if (strcmp(argv[i],"-max_iters") == 0) {
      i++;
      max_iterations = atoi(argv[i]);
    } else if (strcmp(argv[i],"-gmres") == 0) {
      solver = SOLVE_GMRES;
    } else if (strcmp(argv[i],"-fgmres") == 0) {
      solver = SOLVE_FGMRES;
    } else if (strcmp(argv[i],"-local") == 0) {
      solver = SOLVE_FGMRES;
      pc = LOCAL;
    } else if (strcmp(argv[i],"-diagonal") == 0) {
      pc = DIAGONAL;
    } else if (strcmp(argv[i],"-help") == 0) {
      printHelpAndExit();
    } else if (strcmp(argv[i],"-mesh") == 0) {
      i++;
      mesh_name = argv[i];
      mesh = true;
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
      // printHelpAndExit();
    }
  }

  solver_options.max_iters = max_iterations;
  solver_options.restart = max_iterations;
  // opts.sparse_local = true;
  double tic, toc;
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

  if (mesh) {
    printf("reading mesh from: %s\n",mesh_name);
    MeshIO::readMsh<point_type,source_type>(mesh_name, panels); // , panels);
  } else {
    Triangulation::UnitSphere(panels, recursions);
    // initialiseSphere(panels, charges, recursions); //, ProblemOptions());
  }

  // run case solving for Phi (instead of dPhi/dn)
  if (second_kind)
    for (auto& it : panels) it.switch_BC();

  // set constant Phi || dPhi/dn for each panel
  charges.resize(panels.size());
  // set up a more complicated charge, from BEM++
  for (unsigned i=0; i<panels.size(); i++) {
#if BEMCPP_TEST
    auto center = panels[i].center;
    double x = center[0], y = center[1], z = center[2];
    double r = norm(center);
    charges[i] = 2*x*z/(r*r*r*r*r) - y/(r*r*r);
#else
    charges[i] = 1.;
#endif
  }
  // charges = std::vector<charge_type>(panels.size(),1.);

  // Build the FMM structure
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, panels, opts);

  // generate the RHS and initial condition
  std::vector<charge_type> x(panels.size(),0.);

  tic = get_time();
  std::vector<result_type> b(panels.size(),0.);
  double tic2, toc2;
  // generate RHS using temporary FMM plan
  {
    tic2 = get_time();
    for (auto& it : panels) it.switch_BC();
    toc2 = get_time();
    printf("Flipping BC: %g\n",toc2-tic2);
    tic2 = get_time();
    FMM_plan<kernel_type> rhs_plan = FMM_plan<kernel_type>(K,panels,opts);
    toc2 = get_time();
    printf("Creating plan: %g\n",toc2-tic2);
    tic2 = get_time();
    b = rhs_plan.execute(charges);
    toc2 = get_time();
    printf("Executing plan: %g\n",toc2-tic2);
    for (auto& it : panels) it.switch_BC();
  }

  toc = get_time();
  double setup_time = toc-tic;

  // Solve the system using GMRES
  // generate the Preconditioner
  tic = get_time();
  /*
  Preconditioners::Diagonal<charge_type> M(K,
                                           plan.source_begin(),
                                           plan.source_end()
                                          );
  // M.print();
  SolverOptions inner_options(1e-2,1,2);
  inner_options.variable_p = true;
  // Preconditioners::FMGMRES<FMM_plan<kernel_type>,Preconditioners::Diagonal<charge_type>> inner(plan, b, inner_options, M);

  // Local preconditioner
  Preconditioners::LocalInnerSolver<FMM_plan<kernel_type>, Preconditioners::Diagonal<result_type>> local(K, panels, b);

  // block diagonal preconditioner
  Preconditioners::BlockDiagonal<FMM_plan<kernel_type>> block_diag(K,panels);

  */
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
  FGMRESContext<result_type> context(x.size(), solver_options.restart);

  if (second_kind) printf("2nd-kind equation being solved\n");
  else             printf("1st-kind equation being solved\n");
#if 1
  if (solver == SOLVE_GMRES && pc == IDENTITY){
    printf("Solver: GMRES\nPreconditioner: Identity\n");
    // straight GMRES, no preconditioner
    // DirectMV<kernel_type> MV(K, panels, panels);
    GMRES(plan,x,b,solver_options);
  }
#else
  else if (solver == SOLVE_GMRES && pc == DIAGONAL) {
    printf("Solver: GMRES\nPreconditioner: Diagonal\n");
    // GMRES, diagonal preconditioner
    GMRES(plan,x,b,solver_options, M, context);
  }
  else if (solver == SOLVE_FGMRES && pc == IDENTITY) {
    printf("Solver: FMRES\nPreconditioner: Identity\n");
    // GMRES, diagonal preconditioner
    FGMRES(plan,x,b,solver_options);
  }
  else if (solver == SOLVE_FGMRES && pc == DIAGONAL) {
    printf("Solver: FGMRES\nPreconditioner: Block Diagonal\n");
    // GMRES, diagonal preconditioner
    FGMRES(plan,x,b,solver_options, block_diag, context);
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
#endif
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
#if BEMCPP_TEST
  std::vector<result_type> analytical(panels.size());
  for (unsigned i=0; i<panels.size(); i++) {
    auto center = panels[i].center;
    double x = center[0], y = center[1], z = center[2];
    double r = norm(center);
    analytical[i] = -(-6 * x * z / (r*r*r*r*r*r) + 2 * y / (r*r*r*r));
  }

  auto ai = analytical.begin();
  for (auto xi : x) {
    // printf("approx: %.4g, analytical: %.4g\n",xi,*ai);
    e += (xi-*ai)*(xi-*ai);
    e2 += (*ai)*(*ai);
    ++ai;
  }
#else
  double an = 1.;
  for (auto xi : x) { e += (xi-an)*(xi-an); e2 += an*an; }
#endif

#define EXTERNAL_ERROR
#ifdef EXTERNAL_ERROR
  std::vector<target_type> outside_point(1);
  outside_point[0] = target_type(point_type(3.,3.,3.),point_type(3.,3.,3.),point_type(3.,3.,3.));
  outside_point[0].center = point_type(3.,3.,3.);
  std::vector<result_type> outside_result_1(1);
  std::vector<result_type> outside_result_2(1);
  outside_result_1[0] = 0.;
  outside_result_2[0] = 0.;

  // first layer
  Direct::matvec(K, panels.begin(), panels.end(), x.begin(), outside_point.begin(), outside_point.end(), outside_result_2.begin());
  // for (auto& pi : panels) pi.switch_BC();
  for (auto& op : outside_point) op.switch_BC();
  Direct::matvec(K, panels.begin(), panels.end(), charges.begin(), outside_point.begin(), outside_point.end(), outside_result_1.begin());
  double exact = 1. / norm(static_cast<point_type>(outside_point[0])) * 1;
  double outside_result = (outside_result_2[0]-outside_result_1[0])/4/M_PI;
  double outside_error = fabs(outside_result-exact)/fabs(exact);
  printf("external phi: %.5g, exact: %.5g, error: %.4e\n",outside_result,exact, outside_error);
#endif

  printf("error: %.3e\n",sqrt(e/e2));
}

