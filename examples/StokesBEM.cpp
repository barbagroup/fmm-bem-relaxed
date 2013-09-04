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
#include "VertFaceReader.hpp"

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
    auto ci = charges[i]; // *( (it->normal[0]<=0.) ? -1 : 1);//*it->Area;
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
  printf("-rbc <int> : number of recursive subdivisions to create a red blood cell - # panels = 2*4^recursions\n");
  printf("-cells <int> : number of red blood cells to generate\n");
  printf("-fixed_p : Disable relaxation\n");
  printf("-help : print this message\n");
  std::exit(0);
}

template <typename SourceType, typename ChargeType>
void initialiseSphere(std::vector<SourceType>& panels,
                      std::vector<ChargeType>&  charges,
                      unsigned recursions = 4)
{
  (void) charges;
  Triangulation::UnitSphere(panels, recursions);
}

template <typename SourceType, typename ChargeType>
void initialiseRBC(std::vector<SourceType>& panels,
                   std::vector<ChargeType>& charges,
                   unsigned recursions = 4)
{
  (void) charges;
  Mat3<double> rotation = Triangulation::RotationMatrix(0.,0.,0.);
  double s[] = {0.,0.,0.};
  Triangulation::RedBloodCell(panels, recursions, rotation, &s[0]);
}

template <typename SourceType, typename ChargeType>
void initialiseMultipleRBC(std::vector<SourceType>& panels,
                           std::vector<ChargeType>& charges,
                           unsigned recursions = 4,
                           unsigned cells = 2)
{
  (void) charges;
  Triangulation::MultipleRedBloodCell(panels, recursions, cells);
}

int main(int argc, char **argv)
{
  int numPanels= 1000, recursions = 4, p = 8, k = 4, kfine = 19, cells = 1;
  double mu = 1e-3;
  char *mesh_name, *vert_name, *face_name;
  bool mesh = false, seperate_mesh = false, rbc = false;
  unsigned p_min = 5;

  //opts.set_mac_theta(0.5);    // Multipole acceptance criteria
  //opts.set_max_per_box(50);
  FMMOptions opts = get_options(argc, argv);
  opts.sparse_local = true;

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
    } else if (strcmp(argv[i],"-pmin") == 0) {
      i++;
      p_min = (unsigned)atoi(argv[i]);
    } else if (strcmp(argv[i],"-k") == 0) {
      i++;
      k = atoi(argv[i]);
    } else if (strcmp(argv[i],"-mu") == 0) {
      i++;
      mu = atof(argv[i]);
    } else if (strcmp(argv[i],"-recursions") == 0) {
      i++;
      recursions = atoi(argv[i]);
    } else if (strcmp(argv[i],"-rbc") == 0) {
      i++;
      recursions = atoi(argv[i]);
      rbc = true;
    } else if (strcmp(argv[i],"-cells") == 0) {
      i++;
      cells = atoi(argv[i]);
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
    } else if (strcmp(argv[i],"-kfine") == 0) {
      i++;
      kfine = atoi(argv[i]);
    } else if (strcmp(argv[i],"-disable_sparse") == 0) {
      opts.sparse_local = false;
    } else if (strcmp(argv[i],"-vert") == 0) {
      i++;
      vert_name = argv[i];
      mesh = true;
      seperate_mesh = true;
    } else if (strcmp(argv[i],"-face") == 0) {
      i++;
      face_name = argv[i];
      mesh = true;
      seperate_mesh = true;
    }/* else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
      printHelpAndExit();
    } */
  }
  solver_options.max_p = p;
  solver_options.p_min = p_min;
  solver_options.max_iters = 100;
  solver_options.restart = 100;

  double tic, toc;
  // Init the FMM Kernel
  typedef StokesSphericalBEM kernel_type;
  kernel_type K(p,k,mu);

  K.set_Kfine(kfine);

  // useful typedefs
  typedef kernel_type::point_type  point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  // Init points and charges
  std::vector<source_type> panels(numPanels);
  std::vector<charge_type> charges(numPanels);

  // surface mesh -- read / create
  if (mesh) {
    if (seperate_mesh) {
      printf("reading mesh from %s, %s\n",vert_name,face_name);
      MeshIO::ReadVertFace<point_type, source_type>(vert_name, face_name, panels);
    } else {
      printf("reading mesh from %s\n",mesh_name);
      MeshIO::readMsh<point_type, source_type>(mesh_name, panels);
    }
  } else {
    if (rbc) {
      // init RBC
      if (cells <= 1) {
        initialiseRBC(panels, charges, recursions);
      } else {
        initialiseMultipleRBC(panels, charges, recursions, cells);
      }
    } else {
      // init sphere
      initialiseSphere(panels, charges, recursions); //, ProblemOptions());
    }
  }

  // set constant Phi || dPhi/dn for each panel
  charges.resize(panels.size());
  // set ux = 1. -- flow over a sphere
  charges = std::vector<charge_type>(panels.size(),charge_type(1.,0.,0.));

  // generate the RHS and initial condition
  std::vector<charge_type> x(panels.size(),charge_type(1.));

  printf("generating RHS\n");
  // generate RHS using direct calculation
  std::vector<result_type> b(panels.size(),result_type(0.));
  tic = get_time();
  {
    // get RHS for velocity operators
    for (auto& it : panels) it.switch_BC();
    FMM_plan<kernel_type> rhs_plan = FMM_plan<kernel_type>(K,panels,opts);
    b = rhs_plan.execute(charges);
    for (auto& it : panels) it.switch_BC();

    double rhs_error = 0.;
    for (auto& bi : b) {
      rhs_error += fabs(bi[0]-4*M_PI)/4/M_PI;
      bi = result_type(4*M_PI,0.,0.);
    }
    printf("rhs error: %.4e\n",rhs_error);
  }

  printf("done\n");
  toc = get_time();
  double setup_time = toc-tic;

  // Build the FMM structure (after RHS for memory purposes)
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, panels, opts);


#if 0
  // print the RHS and exit
  for (auto bi : b) {
    assert(!isnan(bi[0]) && !isnan(bi[1]) && !isnan(bi[2]));
    std::cout << bi << std::endl;
  }
  std::exit(0);
#endif

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

  double area_sum = 0.;
  // total force in x = sum t^j_x*Area_j
  double t_exact = 1.5*mu*1/1;
  printf("t_exact : %.4g\n",t_exact);
  double e = 0., e2 = 0.;
  for (auto& it : panels) {
    fx += x[i][0]*it.Area;
    fy += x[i][1]*it.Area;
    fz += x[i][2]*it.Area;
    area_sum += it.Area;
    double approx = x[i][0]*it.Area;
    e  += (approx-t_exact)*(approx-t_exact);
    e2 += t_exact*t_exact;
  }
  double analytical_soln = 6*M_PI*mu; // Ux = 1, R = 1
  double analytical_area = 4*M_PI;
  printf("\n\nFx: %.4lg, analytical: %.4lg\n",fx,analytical_soln);
  printf("Fy: %.4g, Fz: %.4g\n",fy,fz);
  double drag_error = fabs(analytical_soln-fx)/fabs(analytical_soln);
  printf("error: %.5e\n",drag_error);
  printf("\n\n");
  printf("\tdrag error per panel : %.5e, average panel area: %.5e\n",drag_error / panels.size(), area_sum / panels.size());
  double area_error = fabs(area_sum - analytical_area) / analytical_area;
  printf("\tArea error : %.5e\n",area_error);

  printf("\n\nPOINTWISE ERRORS\n");
  printf("\terror: %.3e\n",sqrt(e/e2));

  auto& K2 = plan.kernel();
  printf("\nTotals:\n");
  std::cout << "\tStokeslet " << K2.stokeslet_str[0] << ", " << K2.stokeslet_str[1] << ", " << K2.stokeslet_str[2] << ", " << K2.stokeslet_str[3] << std::endl;
  std::cout << "\tStresslet " << K2.stresslet_str[0] << ", " << K2.stresslet_str[1] << ", " << K2.stresslet_str[2] << ", " << K2.stresslet_str[3] << std::endl;

#define EXTERNAL_ERROR
#ifdef EXTERNAL_ERROR
  std::vector<target_type> outside_point(1);
  // point on plane of centerline of sphere
  outside_point[0] = target_type(point_type(3.1,3.,0.),point_type(3.,3.1,0.),point_type(3.,3.,0.));
  outside_point[0].center = point_type(3.,3.,0.);
  std::vector<result_type> outside_result_1(1);
  std::vector<result_type> outside_result_2(1);
  outside_result_1[0] = result_type(0.);
  outside_result_2[0] = result_type(0.);

  // first layer
  Direct::matvec(K, panels.begin(), panels.end(), x.begin(), outside_point.begin(), outside_point.end(), outside_result_2.begin());
  // for (auto& pi : panels) pi.switch_BC();
  for (auto& op : outside_point) op.switch_BC();
  Direct::matvec(K, panels.begin(), panels.end(), charges.begin(), outside_point.begin(), outside_point.end(), outside_result_1.begin());

  double theta = atan2(3.,3.);
  printf("theta: %.4g\n",theta);

  double r = norm(static_cast<point_type>(outside_point[0])) * 1;
  double u_r = cos(theta)*(1+1/2./r/r/r - 3./2./r);
  double u_theta = sin(theta)*(-1. + 1./4./r/r/r + 3./4./r);

  double u_x = cos(theta)*u_r - sin(theta)*u_theta;
  double u_y = sin(theta)*u_r + cos(theta)*u_theta;

  result_type exact(u_x, u_y, 0);
  result_type outside_result = (outside_result_2[0]-outside_result_1[0])/4/M_PI;

  printf("u_x_u: %.4g, u_x_t: %.4g\n",outside_result_2[0][0],outside_result_1[0][0]);
  double x_error = fabs(outside_result[0]-exact[0])/fabs(exact[0]);
  double y_error = fabs(outside_result[1]-exact[1])/fabs(exact[1]);
  printf("u_x: %.4g, %.4g,  u_y: %.4g, %.4g\n",u_x, outside_result[0], u_y, outside_result[1]);
  printf("u_x error: %.4e, u_y error: %.4e\n",x_error, y_error);
#endif
}

