#include "executor/INITM.hpp"
#include "executor/INITL.hpp"

#include "LaplaceCartesian.hpp"
#include "LaplaceSpherical.hpp"

#include "StokesSpherical.hpp"

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
  std::vector<source_type> s(1);
  s[0] = source_type(0,0,0);

  // init charge
  std::vector<charge_type> c(1);
  c[0] = charge_type(1,2,3);

  // init target
  std::vector<target_type> t(1);
  t[0] = target_type(0.9,0,0);
  // target_type t = target_type(0.9,0,0); // 1,1);

  // init results vectors for exact, FMM
  std::vector<result_type> rexact(1);
  // result_type rexact;
  result_type rm2p;
  result_type rfmm;

  // test direct
  // rexact = K(t,s) * c;
  K.P2P(s.begin(),s.end(),c.begin(),t.begin(),t.end(),rexact.begin());
  // rexact = result_type(rexact[0],0.,0.,0.);

  // setup intial multipole expansion
  multipole_type M;
  point_type M_center(0.125,0,0); // 0.125,0.125);
  INITM::eval(K, M, M_center, 1u);
  K.P2M(s[0], c[0], M_center, M);

  // test M2P
  K.M2P(M, M_center, t[0], rm2p);

  // test M2L, L2P
  local_type L;
  point_type L_center(0.875,0,0); // 0.875,0.875);
  auto d = L_center - M_center;
  printf("DIST: (%lg, %lg, %lg) : %lg\n",d[0],d[1],d[2],norm(d));
  // L_center = point_type(t);
  INITL::eval(K, L, L_center, 1u);
  K.M2L(M, L, L_center - M_center);
  K.L2P(L, L_center, t[0], rfmm);

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

  std::cout << "FMM L1 rel error: "
            << std::scientific << l1_rel_error(rfmm, rexact) << std::endl;
  std::cout << "FMM L2 error:     "
            << std::scientific << l2_error(rfmm, rexact) << std::endl;
  std::cout << "FMM L2 rel error: "
            << std::scientific << l2_rel_error(rfmm, rexact) << std::endl;
  */
}

int main(int argc, char** argv)
{
  (void) argc;
  (void) argv;

  //LaplaceCartesian<5> K;
  // LaplaceSpherical K(5);
  StokesSpherical K(5);

  single_level_test(K);
}

