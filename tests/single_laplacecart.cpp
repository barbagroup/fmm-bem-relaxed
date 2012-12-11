#include <LaplaceCartesian.hpp>
#include <Direct.hpp>

#include <iostream>
#include <cstring>

template <typename Kernel>
void do_test(const Kernel& K) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::point_type point_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;
  typedef typename kernel_type::multipole_type multipole_type;
  typedef typename kernel_type::local_type local_type;

  // init source
  source_type s = source_type(0,0,0);

  // init charge
  charge_type c = 1;

  // init target
  target_type t = target_type(1,0,0);

  // init results vectors for exact, FMM
  result_type rexact;
  result_type rm2p;
  result_type rfmm;

  // test direct
  rexact = K(t,s) * c;

  // setup intial multipole expansion
  multipole_type M;
  point_type source_center(0.125,0.125,0.125);
  K.P2M(s, c, source_center, M);

  // test M2P
  K.M2P(M, source_center, t, rm2p);

  // test M2L, L2P
  local_type L;
  point_type target_center(0.875,0.125,0.125);
  K.M2L(M, L, target_center - source_center);
  K.L2P(L, target_center, t, rfmm);

  // check errors
  std::cout << "rexact = " << rexact << std::endl;
  std::cout << "rm2p = " << rm2p << std::endl;
  std::cout << "rfmm = " << rfmm << std::endl;

  double M2P_error = fabs(rexact[0] - rm2p[0]) / rexact[0];
  double FMM_error = fabs(rexact[0] - rfmm[0]) / rexact[0];

  printf("M2P error: %.4le\n", M2P_error);
  printf("FMM error: %.4le\n", FMM_error);
}



int main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  do_test(LaplaceCartesian<5>());
}

