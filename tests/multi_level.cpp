#include "executor/INITM.hpp"
#include "executor/INITL.hpp"

#include "LaplaceCartesian.hpp"
#include "LaplaceSpherical.hpp"
#include "YukawaCartesian.hpp"
#include "StokesSpherical.hpp"

#include "Math.hpp"

#include <iostream>

// #define TREECODE_ONLY

template <typename Kernel>
void two_level_test(const Kernel& K)
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
  std::vector<source_type> s(1);
  s[0] = source_type(0,0,0);

  // init charge
  std::vector<charge_type> c(1);
  c[0] = charge_type(1,2,3);

  // init target
  std::vector<target_type> t(1);
  t[0] = target_type(0.98,0.98,0.98);

  // init results vectors for exact, FMM
  std::vector<result_type> rexact(1);
  result_type rm2p;
  result_type rfmm;

  // test direct
  K.P2P(s.begin(),s.end(),c.begin(),t.begin(),t.end(),rexact.begin());

  // setup intial multipole expansion
  multipole_type M;
  const point_type M_center(0.05, 0.05, 0.05);
  INITM::eval(K, M, M_center, 2u);
  K.P2M(s[0], c[0], M_center, M);

  // perform M2M
  multipole_type M2;
  const point_type M2_center(0.1, 0.1, 0.1);
  INITM::eval(K, M2, M2_center, 1u);
  K.M2M(M, M2, M2_center - M_center);
  //K.P2M(s,c,M2_center,M2);

  // test M2P
  K.M2P(M2, M2_center, t[0], rm2p);

  // test M2L, L2P
#ifndef TREECODE_ONLY
  local_type L2;
  point_type L2_center(0.9, 0.9, 0.9);
  INITL::eval(K, L2, L2_center, 1u);
  K.M2L(M2, L2, L2_center - M2_center);

  // test L2L
  local_type L;
  point_type L_center(0.95, 0.95, 0.95);
  INITL::eval(K, L, L_center, 2u);
  K.L2L(L2, L, L_center - L2_center);

  // test L2P
  K.L2P(L2, L2_center, t[0], rfmm);
#endif

  // check errors
  std::cout << "rexact = " << rexact[0] << std::endl;
  std::cout << "rm2p = " << rm2p << std::endl;
  std::cout << "rfmm = " << rfmm << std::endl;

/*
  std::cout << "M2P L1 rel error: "
            << std::scientific << l1_rel_error(rm2p, rexact) << std::endl;
  std::cout << "M2P L2 error:     "
            << std::scientific << l2_error(rm2p, rexact) << std::endl;
  std::cout << "M2P L2 rel error: "
            << std::scientific << l2_rel_error(rm2p, rexact) << std::endl;
#ifndef TREECODE_ONLY
  std::cout << "FMM L1 rel error: "
            << std::scientific << l1_rel_error(rfmm, rexact) << std::endl;
  std::cout << "FMM L2 error:     "
            << std::scientific << l2_error(rfmm, rexact) << std::endl;
  std::cout << "FMM L2 rel error: "
            << std::scientific << l2_rel_error(rfmm, rexact) << std::endl;
#endif
*/
}


int main(int argc, char **argv)
{
  (void) argc;
  (void) argv;

  //LaplaceCartesian<5> K;
  //LaplaceSpherical K(10);
  //YukawaCartesian K(10, 0.1);
  StokesSpherical K(5);

  two_level_test(K);
}

