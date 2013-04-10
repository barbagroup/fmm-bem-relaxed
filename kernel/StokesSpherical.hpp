#pragma once

/** @file StokesSpherical.hpp
 * @brief Implements Stokes kernel with repeated Laplace calls
 */

#include "LaplaceSpherical.hpp"

class StokesSpherical : public LaplaceSpherical
{
 public:
  //! point type
  typedef LaplaceSpherical::point_type point_type;
  //! source type
  typedef LaplaceSpherical::source_type source_type;
  //! target type
  typedef LaplaceSpherical::target_type target_type;
  //! charge type - { f1, f2, f3 }
  typedef Vec<3, LaplaceSpherical::charge_type> charge_type;
  //! return of kernel evaluation
  typedef LaplaceSpherical::kernel_value_type kernel_value_type;
  //! product of charge and KVT
  typedef Vec<3, real> result_type;

  //! multipole type
  typedef Vec<4, LaplaceSpherical::multipole_type> multipole_type;
  //! local type
  typedef Vec<4, LaplaceSpherical::local_type> local_type;

  //! default (delegating) constructor
  StokesSpherical() : StokesSpherical(5) {};
  //! Constructor
  StokesSpherical(int p) : LaplaceSpherical(p) {};

  /** Initialize a multipole expansion */
  void init_multipole(multipole_type& M, const point_type extents, unsigned level) const {
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::init_multipole(M[i], extents, level);
    }
  }
  /** Initialize a local expansion */
  void init_local(local_type& L, const point_type extents, unsigned level) const {
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::init_local(M[i], extents, level);
    }
  }

  /** Kernel evaluation
   * K(t,s)
   * this is kind of useless -- I need to define a P2P anyway, as charge * KVT doesn't work for this app...
   */
  kernel_value_type operator()(const target_type& t, const source_type& s) const {
    // do shit here
    point_type dist = s - t;
    real R2 = norm(dist);
    real invR2 = 1.0 / R2;
    if (R2 < 1e-8) invR2 = 0;
    real invR = std::sqrt(invR2)
  }

  template <typename SourceIter, typename ChargeIter,
            typename TargetIter, typename ResultIter>
  void P2P(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first) const {
    // do crap here
    for ( ; t_first != t_last; ++t_first, ++r_first) {
      for ( ; s_first != s_last; ++s_first, ++c_first) {
        // interact target, source & charge here
        point_type dist = *t_first - *s_first;
        auto r2 = norm(dist);
        auto invR2 = 1./r2;
        if (invR2 < 1e-8) invR2 = 0;
        auto invR = std::sqrt(invR2);
        auto H = invR * invR2; // 1 / R^3

        auto c = *c_first;
        auto fdx = dist[0]*c[0] + dist[1]*c[1] + dist[2]*c[2];

        auto r = *r_first;

        r[0] += H*(c[0]*r2 + fdx*dist[0];
        r[1] += H*(c[1]*r2 + fdx*dist[1];
        r[2] += H*(c[2]*r2 + fdx*dist[2];
      }
    }
  }

  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    (void) source;
    (void) charge;
    (void) center;
    (void) M;
    // modifications needed here
  }

  void M2M(const multipole_type& Msource,
                 multipole_type& Mtarget, const point_type& translation) const {
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::M2M(Msource[i], Mtarget[i], translation);
    }
  }

  void M2P(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result) const {
    //
    (void) M;
    (void) center;
    (void) target;
    (void) result;
  }

  void M2L(const multipole_type& Msource,
                 local_type& Ltarget, const point_type& translation) const {
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::M2L(Msource[i], Ltarget[i], translation);
    }
  }

  void L2L(const local_type& Lsource,
                 local_type& Ltarget, const point_type& translation) const {
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::M2L(Lsource[i], Ltarget[i], translation);
    }
  }

  void L2P(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
    (void) L;
    (void) center;
    (void) target;
    (void) result;
  }

}; // end class StokesSpherical
