/** @file serialBEM.cpp
 * @brief Testing and debugging script for FMM-BEM applications
 */

#include <fstream>
// general includes
#include "FMM_plan.hpp"
#include "DirectMatvec.hpp"
#include "Triangulation.hpp"
#include "SolverOptions.hpp"

// kernel
#include "StokesSphericalBEM.hpp"

// solvers and preconditioners
#include "GMRES_Stokes.hpp"
#include "LocalPC_Stokes.hpp"
#include "BlockDiagonalPC_Stokes.hpp"

// mesh reader
#include "MshReader.hpp"

// timing
#include "timing.hpp"

enum SOLVERS {
               SOLVE_GMRES,
               SOLVE_FGMRES
             };
enum PRECONDITIONERS {
                       IDENTITY,
                       DIAGONAL, // block-diagonal solve
                       LOCAL     // full local solve
                     };

struct ProblemOptions
{
  typedef enum { VELOCITY_SET, TRACTION_SET } BoundaryCondition;
  BoundaryCondition bc_ = VELOCITY_SET;
  double value_ = 1.;
  int recursions = 4;

  ProblemOptions() : bc_(VELOCITY_SET), value_(1.) {};
  ProblemOptions(double value) : bc_(VELOCITY_SET), value_(value) {};
  ProblemOptions(BoundaryCondition bc, double value) : bc_(bc), value_(value) {};

  double getValue() { return value_; };
  BoundaryCondition getBC() { return bc_; };
};

template <typename Panel, typename Charge>
void output_solution(const std::vector<Panel>& panels, const std::vector<Charge>& charges)
{
  // output the vertices, faces and value
  std::ofstream face("out.face");
  std::ofstream vert("out.vert");
  std::ofstream c("out.charge");

  int i=0;
  int vnum = 1;
  for (auto it=panels.begin(); it!=panels.end(); ++it) {
    vert << it->vertices[0][0] << "," << it->vertices[0][1] << "," << it->vertices[0][2] << std::endl;
    vert << it->vertices[1][0] << "," << it->vertices[1][1] << "," << it->vertices[1][2] << std::endl;
    vert << it->vertices[2][0] << "," << it->vertices[2][1] << "," << it->vertices[2][2] << std::endl;
    face << vnum << "," << vnum+1 << "," << vnum+2 << std::endl;
    auto ci = charges[i]; // *it->Area;
    //c << ci[0]*it->normal[0] << "," << ci[1]*it->normal[0] << "," << ci[2]*it->normal[0] << std::endl;
    c << ci[0] << "," << ci[1] << "," << ci[2] << std::endl;
    vnum+=3;
    i++;
  }
  face.close();
  vert.close();
  c.close();
}

void printHelpAndExit()
{
  printf("StokesBEM : FMM-BEM for Stokes problems\n");
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
  int numPanels= 1000, recursions = 4, p = 8, k = 4;
  double mu = 1e-3;
  char *mesh_name;
  bool mesh = false;

  FMMOptions opts = get_options(argc, argv);
  opts.set_mac_theta(0.5);    // Multipole acceptance criteria
  opts.set_max_per_box(50);

  SolverOptions solver_options;
  SOLVERS solver = SOLVE_GMRES;
  PRECONDITIONERS pc = IDENTITY;

  // parse command line args
  // check if no arguments given
  if (argc == 1) printHelpAndExit();
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-p") == 0) {
      i++;
      p = atoi(argv[i]);
    } else if (strcmp(argv[i],"-k") == 0) {
      i++;
      k = atoi(argv[i]);
    } else if (strcmp(argv[i],"-recursions") == 0) {
      i++;
      recursions = atoi(argv[i]);
    } else if (strcmp(argv[i],"-fixed_p") == 0) {
      solver_options.variable_p = false;
    } else if (strcmp(argv[i],"-solver_tol") == 0) {
      i++;
      solver_options.residual = atof(argv[i]);
    } else if (strcmp(argv[i],"-mesh") == 0) {
      i++;
      mesh_name = argv[i];
      mesh = true;
    } else if (strcmp(argv[i],"-fgmres") == 0) {
      solver = SOLVE_FGMRES;
    } else if (strcmp(argv[i],"-diagonal") == 0 || strcmp(argv[i],"-diag") == 0) {
      solver = SOLVE_FGMRES;
      pc = DIAGONAL;
    } else if (strcmp(argv[i],"-local") == 0) {
      solver = SOLVE_FGMRES;
      pc = LOCAL;
    } else if (strcmp(argv[i],"-help") == 0) {
      printHelpAndExit();
    }/* else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
      printHelpAndExit();
    } */
  }

  double tic, toc;
  tic = get_time();
  // Init the FMM Kernel
  typedef StokesSphericalBEM kernel_type;
  kernel_type K(p,k,mu);

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
    printf("reading mesh from %s\n",mesh_name);
    readMsh<point_type, source_type>(mesh_name, panels);
  } else {
    initialiseSphere(panels, charges, recursions); //, ProblemOptions());
  }

  // set constant Phi || dPhi/dn for each panel
  charges.resize(panels.size());
  // set ux = 1. -- flow over a sphere
  charges = std::vector<charge_type>(panels.size(),charge_type(1.,0.,0.));

  // Build the FMM structure
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, panels, opts);

  // generate the RHS and initial condition
  std::vector<charge_type> x(panels.size(),charge_type(0.));

  // generate RHS using direct calculation
  std::vector<result_type> b(panels.size(),result_type(0.));
  {
    // get RHS for velocity operators
    for (auto& it : panels) it.switch_BC();
    FMM_plan<kernel_type> rhs_plan = FMM_plan<kernel_type>(K,panels,opts);
    b = rhs_plan.execute(charges);
    for (auto& it : panels) it.switch_BC();
  }
  toc = get_time();
  double setup_time = toc-tic;

  // print the RHS and exit
  /*
  for (auto bi : b) {
    assert(!isnan(bi[0]) && !isnan(bi[1]) && !isnan(bi[2]));
    std::cout << bi << std::endl;
  }
  std::exit(0);
  */

  // Solve the system using GMRES
  // generate the Preconditioner
  //Preconditioners::Diagonal<charge_type> M(K,panels.begin(),panels.end());
  // FGMRES(plan, x, b, SolverOptions());
  //fmm_gmres(plan, x, b, SolverOptions());
  //direct_gmres(K, panels, x, b, SolverOptions());

  tic = get_time();
  if (solver == SOLVE_GMRES) {
    printf("Solver: GMRES, Preconditioner: Identity\n");
    GMRES(plan, x, b, solver_options);
  } else if (solver == SOLVE_FGMRES) {
    if (pc == IDENTITY) {
      printf("Solver: FGMRES, Preconditioner: Identity\n");
      FGMRES(plan, x, b, solver_options);
    } else if (pc == DIAGONAL) {
      printf("Solver: FGMRES, Preconditioner: Block-Diagonal\n");
      Preconditioners::BlockDiagonal<FMM_plan<kernel_type>> block_diag(K, panels);
      FGMRES(plan, x, b, solver_options, block_diag);
    } else if (pc == LOCAL) {
      printf("Solver: FGMRES, Preconditioner: Local Solve\n");
      Preconditioners::LocalInnerSolver<FMM_plan<kernel_type>> local(K, panels, b);
      FGMRES(plan, x, b, solver_options, local);
    } else {
      printf("No known preconditioner specified\n");
      std::exit(0);
    }
  } else {
    printf("No known solver specified\n");
    std::exit(0);
  }
  toc = get_time();
  double solve_time = toc-tic;

  printf("\nTIMING:\n");
  printf("\tsetup : %.4es\n",setup_time);
  printf("\tsolve : %.4es\n\n",solve_time);

  double fx = 0., fy = 0., fz = 0.;
  int i=0;

  // output the solution
  output_solution(panels, x);

  // total force in x = sum t^j_x*Area_j
  for (auto& it : panels) {
    fx += x[i][0]*it.Area;
    fy += x[i][1]*it.Area;
    fz += x[i][2]*it.Area;
  }
  double analytical_soln = -6*M_PI*mu; // Ux = 1, R = 1
  printf("\n\nFx: %.4lg, analytical: %.4lg\n",fx,analytical_soln);
  printf("Fy: %.4g, Fz: %.4g\n",fy,fz);
  printf("error: %.5e\n",fabs(analytical_soln-fx)/fabs(analytical_soln));
}

