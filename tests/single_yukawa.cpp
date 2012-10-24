// #include <SphericalLaplaceKernel.hpp>
#include <CartesianYukawaKernel.hpp>
#include <Direct.hpp>
#include <cstring>

int main(int argc, char **argv)
{
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

  // init source
  std::vector<point_type> points(1);
  points[0] = point_type(0,0,0);
  //points[1] = point_type(0.05,0,0);

  // init charge 
  std::vector<charge_type> charges(1);
  charges[0] = 1.;
  //charges[1] = 1.;

  // init target
  std::vector<point_type> target(1);
  target[0] = point_type(1,1,1);
  //target[1] = point_type(1,1,0.95);
  //target[2] = point_type(1,0.95,1);
  //target[3] = point_type(0.95,1,1);

  // init results vectors for exact, FMM
  std::vector<result_type> exact(target.size());
  std::vector<result_type> M2P(target.size());

  // init source, target M, L arrays
  multipole_type M1, M2;
  K.init_multipole(M1,0);
  K.init_multipole(M2,0);

  // setup intial multipole expansion
  point_type source_center(0.05,0.05,0.05);
  K.P2M(points.begin(),points.end(),charges.begin(),source_center,M1);
  // second Multipole expansion
  // K.P2M(points.begin()+1, points.end(),charges.begin()+1,source_center,M2);

  // test direct
  Direct::matvec(K,points.begin(),points.end(),charges.begin(),target.begin(),target.end(),exact.begin());

  // test M2P
  K.M2P(M1,source_center,target.begin(),target.end(),M2P.begin());
  // second eval
  // K.M2P(M2,source_center,target.begin(),target.end(),M2P.begin());

  // check errors
  double M2P_error  = fabs(exact[0][1] - M2P[0][1]) / exact[0][1];
         //M2P_error += fabs(exact[1][1] - M2P[1][1]) / exact[1][1];
         //M2P_error += fabs(exact[2][1] - M2P[2][1]) / exact[2][1];
         //M2P_error += fabs(exact[3][1] - M2P[3][1]) / exact[3][1];

  printf("exact: %.3lg, M2P: %.3lg\n",exact[0][0],M2P[0][0]);
  printf("exact: %.3lg, M2P: %.3lg\n",exact[0][1],M2P[0][1]);
  printf("exact: %.3lg, M2P: %.3lg\n",exact[0][2],M2P[0][2]);
  printf("exact: %.3lg, M2P: %.3lg\n",exact[0][3],M2P[0][3]);
  printf("M2P error: %.4lg\n",M2P_error);

}

