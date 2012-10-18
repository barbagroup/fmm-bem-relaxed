#include <SphericalLaplaceKernel.hpp>
#include <Direct.hpp>
#include <cstring>

int main(int argc, char **argv)
{
  typedef SphericalLaplaceKernel kernel_type;
  typedef kernel_type::point_type point_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;
  typedef kernel_type::multipole_type multipole_type;
  typedef kernel_type::local_type local_type;

  int P=5;

  for (int i=1; i<argc; i++)
  {
    if (strcmp(argv[i],"-P")==0) {
      i++;
      P = atoi(argv[i]);
    }
    else {
      printf("[W]: unknown command line argument: \"%s\"\n",argv[i]);
    }
  }
  kernel_type K(P);

  // init source
  std::vector<point_type> points(1);
  points[0] = point_type(0,0,0);

  // init charge 
  std::vector<charge_type> charges(1);
  charges[0] = 1.;

  // init target
  std::vector<point_type> target(1);
  target[0] = point_type(1,1,1);

  // init results vectors for exact, FMM
  std::vector<result_type> exact(target.size());
  std::vector<result_type> M2P(target.size());
  std::vector<result_type> FMM(target.size());

  // init source, target M, L arrays
  multipole_type M;
  K.init_multipole(M,0);
  local_type L;
  K.init_local(L,0);

  // setup intial multipole expansion
  point_type source_center(0.05,0.05,0.05);
  K.P2M(points.begin(),points.end(),charges.begin(),source_center,M);

  // target 'box' center
  point_type target_center(0.95,0.95,0.95);

  // test direct
  Direct::matvec(K,points.begin(),points.end(),charges.begin(),target.begin(),target.end(),exact.begin());

  // test M2P
  K.M2P(M,source_center,target.begin(),target.end(),M2P.begin());

  // test M2L, L2P
  auto translation = target_center - source_center;
  K.M2L(M,L,translation);

  K.L2P(L,target_center,target.begin(),target.end(),FMM.begin());

  // check errors
  double M2P_error = fabs(exact[0][0] - M2P[0][0]) / exact[0][0];
  double FMM_error = fabs(exact[0][0] - FMM[0][0]) / exact[0][0];

  printf("M2P error: %.4lg\n",M2P_error);
  printf("FMM error: %.4lg\n",FMM_error);


}

