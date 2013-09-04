#include <FMM_plan.hpp>
#include <LaplaceSpherical.hpp>

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
  opts.set_max_per_box(126);

  std::vector<std::pair<unsigned,double>> times;
  double tic, toc;

  for (double n = 4; n <= 6; n += 0.125) {
    int numBodies = int(pow(10,n));

    std::vector<point_type> points(numBodies);
    for (int k=0; k<numBodies; ++k) {
      points[k] = point_type(drand(),drand(),drand());
    }

    std::vector<charge_type> charges(numBodies);
    for (int k=0; k<numBodies; ++k) {
      charges[k] = drand();
    }

    FMM_plan<kernel_type> plan = FMM_plan<kernel_type>(K, points, opts);

    // initialise plan & solve
    tic = get_time();
    auto result = plan.execute(charges);
    toc = get_time();

    printf("%d\t%.3e\n", numBodies, toc-tic);
    times.push_back(std::make_pair(numBodies,toc-tic));
  }

  for (auto it=times.begin(); it!=times.end(); ++it) {
    printf("%d\t%.3e\n", it->first, it->second);
  }

  for (auto it=times.begin(); it!=times.end(); ++it) {
    printf("%d\n",it->first);
  }
  for (auto it=times.begin(); it!=times.end(); ++it) {
    printf("%.3e\n",it->second);
  }
}

