#include <FMM_plan.hpp>
#include <LaplaceSpherical.hpp>
#include <cmath>
#include <numeric>

inline double drand()
{
  return ::drand48();
}

int main()
{
  typedef LaplaceSpherical kernel_type;
  kernel_type K(5);
  typedef kernel_type::point_type point_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  FMMOptions opts;
  opts.set_mac_theta(.5);
  opts.set_max_per_box(125);    // optimal ncrit

  std::vector<std::pair<unsigned,double>> times;
  double tic, toc;

  
  int numBodies = 10000;
  // initialize points
  std::vector<point_type> points(numBodies);
  for (int k=0; k<numBodies; ++k){
          points[k] = point_type(drand(), drand(), drand());
  }

  // initialize charges
  std::vector<charge_type> charges(numBodies);
  for (int k=0; k<numBodies; ++k){
          charges[k] = drand();
  }

  // create FMM plan
  FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, points, opts);
  // execute FMM
  // run 3 times and make an average
  int nt = 3;    // number of identical runs for timing
  std::vector<double> timings(nt);
  std::vector<result_type> result(numBodies);
  for (int i=0; i<nt; i++){
	tic = get_time();
  	result = plan.execute(charges);
  	toc = get_time();
	timings[i] = toc-tic;
  }
  
  double FMM_time = std::accumulate(timings.begin(), timings.end(), 0.0) / timings.size();
  std::cout << "FMM execution time: " << FMM_time << std::endl;
  // init exact
  std::vector<result_type> exact(numBodies);
  // compute using direct summation
  tic = get_time();
  Direct::matvec(K, points, charges, exact);
  toc = get_time();
 
  std::cout << "direct summation execution time: " << toc-tic << std::endl;
  
  // calculate l2 norm of relative error of the Laplace force
  double e1 = 0;    // numerator
  double e2 = 0;    // denominator
  for (int k=0; k<numBodies; ++k){
	for (int m=1; m<4; ++m){
		e1 += (result[k][m] - exact[k][m]) * (result[k][m] - exact[k][m]);
		e2 += exact[k][m] * exact[k][m];
	}
  }
  std::cout << "the l2-norm of the relative error of Laplace force: " << sqrt(e1/e2) << std::endl;
  return 0;

}
