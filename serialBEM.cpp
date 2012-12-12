/** @file serialBEM.cpp
 * @brief Testing and debugging script for FMM-BEM applications
 */


#include "FMM_plan.hpp"
#include "LaplaceSphericalBEM.hpp"
#include "Triangulation.hpp"
#include "gmres.hpp"

struct SolverOptions
{
  double residual;
  int max_iters, restart;

  SolverOptions() : residual(1e-5), max_iters(50), restart(50) {};
};

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
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
    }
  }

  // Init the FMM Kernel
  typedef LaplaceSphericalBEM kernel_type;
  kernel_type K(p,k);

  typedef kernel_type::point_type  point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  
  // Init points and charges
  std::vector<source_type> panels(numPanels);
  std::vector<charge_type> charges(numPanels);
  initialiseSphere(panels, charges, recursions);

  // run case solving for Phi (instead of dPhi/dn)
  for (auto& it : panels) it.switch_BC();
  charges.resize(panels.size());
  charges = std::vector<charge_type>(panels.size(),1.);

  // Build the FMM structure
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, panels, opts);

  // generate the RHS and initial condition
  std::vector<charge_type> x(panels.size(),0.);

  // generate RHS using direct calculation
  for (auto& it : panels) it.switch_BC();
  std::vector<result_type> b(panels.size(),0.);
  std::vector<result_type> b2(panels.size(),0.);
  b2 = plan.execute(charges);
  Direct::matvec(K,panels,charges,b);
  // b = plan.execute(charges);
  for (auto& it : panels) it.switch_BC();

  // Solve the system using GMRES
  fmm_gmres(plan, x, b, SolverOptions());
  //direct_gmres(K, panels, x, b, SolverOptions());

  // check errors -- analytical solution for dPhi/dn = 1.
  double e = 0.;
  double e2 = 0.;
  double an = 1.;
  for (auto xi : x) { e += (xi-an)*(xi-an); e2 += an*an; }
  
  printf("error: %.3e\n",sqrt(e/e2));

}

