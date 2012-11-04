#include <CartesianYukawaKernel.hpp>
#include <UnitKernel.hpp>

#include <Direct.hpp>
#include <cstring>

int main(int argc, char **argv)
{
  // typedefs
  typedef CartesianYukawaKernel kernel_type;
  typedef kernel_type::point_type point_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;
  typedef kernel_type::multipole_type multipole_type;
  typedef kernel_type::local_type local_type;

  int P=2;

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
  kernel_type K(P,0.);
  //kernel_type K;

  // init source
  std::vector<point_type> points(1);
  points[0] = point_type(0.12,0.12,0.12);

  // init charge
  std::vector<charge_type> charges(1);
  charges[0] = 1.;

  // init target
  std::vector<point_type> target(1);
  target[0] = point_type(0.9,0.9,0.9);

  // init results vectors for exact, FMM
  std::vector<result_type> exact(target.size());
  std::vector<result_type> M2P(target.size());

  // init source, target M, L arrays
  multipole_type M, M2;
  K.init_multipole(M,0);
  K.init_multipole(M2,0);

  // setup intial multipole expansion
  point_type source_center(0.125,0.125,0.125);
  K.P2M(points.begin(),points.end(),charges.begin(),source_center,M);

  // perform M2M
  point_type M_trans_center(0.25,0.25,0.25);
  auto M2M_translation = M_trans_center - source_center;
  K.M2M(M,M2,M2M_translation);

  // test direct
  Direct::matvec(K,points.begin(),points.end(),charges.begin(),target.begin(),target.end(),exact.begin());

  // test M2P
  K.M2P(M2,M_trans_center,target.begin(),target.end(),M2P.begin());

  // check errors
  double M2P_error = fabs(exact[0][0] - M2P[0][0]) / exact[0][0];

  for (unsigned i=0; i<M2.size(); i++) {
    printf("M[%02d]: %lg\n",i,M2[i]);
  }

  printf("exact: %.3lg, M2P: %.3lg\n",exact[0][0],M2P[0][0]);
  printf("M2P error: %.4lg\n",M2P_error);
}

