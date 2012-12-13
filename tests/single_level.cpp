#include "INITM.hpp"
#include "INITL.hpp"

#include "LaplaceCartesian.hpp"
#include "LaplaceSpherical.hpp"
#include "YukawaCartesian.hpp"

#include "Math.hpp"

#include <iostream>

template <typename Kernel>
void single_level_test(const Kernel& K)
{
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
  point_type M_center(0.125,0.125,0.125);
  INITM::eval(K, M, 0.25);
  K.P2M(s, c, M_center, M);

  // test M2P
  K.M2P(M, M_center, t, rm2p);

  // test M2L, L2P
  local_type L;
  INITL::eval(K, L, 0.25);
  point_type L_center(0.875,0.125,0.125);
  K.M2L(M, L, L_center - M_center);
  K.L2P(L, L_center, t, rfmm);

  // check errors
  std::cout << "rexact = " << rexact << std::endl;
  std::cout << "rm2p = " << rm2p << std::endl;
  std::cout << "rfmm = " << rfmm << std::endl;

  std::cout << "M2P L1 rel error: "
            << std::scientific << l1_rel_error(rm2p, rexact) << std::endl;
  std::cout << "M2P L2 error:     "
            << std::scientific << l2_error(rm2p, rexact) << std::endl;
  std::cout << "M2P L2 rel error: "
            << std::scientific << l2_rel_error(rm2p, rexact) << std::endl;

  std::cout << "FMM L1 rel error: "
            << std::scientific << l1_rel_error(rfmm, rexact) << std::endl;
  std::cout << "FMM L2 error:     "
            << std::scientific << l2_error(rfmm, rexact) << std::endl;
  std::cout << "FMM L2 rel error: "
            << std::scientific << l2_rel_error(rfmm, rexact) << std::endl;
}

int main(int argc, char** argv)
{
  (void) argc;
  (void) argv;

  //LaplaceCartesian<5> K;
  LaplaceSpherical K(5);
  //YukawaCartesian K(5, 0);

  single_level_test(K);
}

