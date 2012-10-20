#include <FMM_plan.hpp>
#include <SphericalLaplaceKernel.hpp>

double get_time()
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (double)(tv.tv_sec+tv.tv_usec*1e-6);
}

inline double drand()
{
  return ::drand48();
}

int main(int argc, char **argv)
{
  typedef SphericalLaplaceKernel kernel_type;
  kernel_type K(5);
  typedef kernel_type::point_type point_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  FMMOptions opts;
  opts.set_theta(.5);
  opts.NCRIT = 150;

  std::vector<std::pair<unsigned,double>> times;
  double tic, toc;

  for (unsigned it=0; it!=25; ++it)
  {
    int numBodies = int(pow(10,(it+24)/8.));

    std::vector<point_type> points(numBodies);
    for (int k=0; k<numBodies; ++k) {
      points[k] = point_type(drand(),drand(),drand());
    }
    auto jpoints = points;

    std::vector<charge_type> charges(numBodies);
    for (int k=0; k<numBodies; ++k) {
      charges[k] = drand();
    }

    // initialise plan & solve
    tic = get_time();
    FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, points, opts);
    auto result = plan.execute(charges,jpoints);
    toc = get_time();
    
    times.push_back(std::make_pair(numBodies,toc-tic));
  }

  for (auto it=times.begin(); it!=times.end(); ++it) {
    printf("%d : %.3e\n",it->first,it->second);
  }
}

