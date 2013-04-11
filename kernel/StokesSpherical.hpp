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
  typedef Vec<3, real> kernel_value_type;
  //! product of charge and KVT
  typedef Vec<3, real> result_type;

  //! multipole type
  //typedef Vec<4, LaplaceSpherical::multipole_type> multipole_type;
  typedef std::vector<LaplaceSpherical::multipole_type> multipole_type;
  //! local type
  //typedef Vec<4, LaplaceSpherical::local_type> local_type;
  typedef std::vector<LaplaceSpherical::local_type> local_type;

  //! default (delegating) constructor
  StokesSpherical() : StokesSpherical(5) {};
  //! Constructor
  StokesSpherical(int p) : LaplaceSpherical(p) {};

  /** Initialize a multipole expansion */
  void init_multipole(multipole_type& M, const point_type& extents, unsigned level) const {
    M = std::vector<LaplaceSpherical::multipole_type>(4);
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::init_multipole(M[i], extents, level);
    }
  }
  /** Initialize a local expansion */
  void init_local(local_type& L, const point_type& extents, unsigned level) const {
    L.resize(4);
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::init_local(L[i], extents, level);
    }
  }

  /** Kernel evaluation
   * K(t,s)
   * this is kind of useless -- I need to define a P2P anyway, as charge * KVT doesn't work for this app...
   */
  /*
  kernel_value_type operator()(const target_type& t, const source_type& s) const {
    // do shit here
    (void) t;
    (void) s;
    //oint_type dist = s - t;
    //real R2 = norm(dist);
    //real invR2 = 1.0 / R2;
    //if (R2 < 1e-8) invR2 = 0;
    // real invR = std::sqrt(invR2);
    return kernel_value_type(0,0,0);
  }
  */

  template <typename SourceIter, typename ChargeIter,
            typename TargetIter, typename ResultIter>
  void P2P(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first) const {
    // do crap here
    auto ti = t_first;
    auto ri = r_first;

    for ( ; ti != t_last; ++ti, ++ri) {
      auto si = s_first;
      auto ci = c_first;
      for ( ; si != s_last; ++si, ++ci) {
        // interact target, source & charge here
        point_type dist = *ti - *si;
        auto r2 = normSq(dist);
        real R1 = r2;
        real R2 = R1;
        real invR = 1. / R1;
        if (r2 < 1e-8) invR = 0;

        //auto invR2 = 1./r2;
        //if (invR2 < 1e-8) invR2 = 0;
        //auto invR = std::sqrt(invR2);
        auto H = std::sqrt(invR) * invR; // 1 / R^3

        //auto c = *c_first;
        auto fdx = dist[0]*(*ci)[0] + dist[1]*(*ci)[1] + dist[2]*(*ci)[2];

        // auto& r = *r_first;

        //r[0] += 1.; // H * (c[0] * R2 + fdx * dist[0]);
        //r[1] += 1.; // H * (c[1] * R2 + fdx * dist[1]);
        //r[2] += 1.; // H * (c[2] * R2 + fdx * dist[2]);
        (*ri)[0] +=  H * ((*ci)[0] * R2 + fdx * dist[0]);
        (*ri)[1] +=  H * ((*ci)[1] * R2 + fdx * dist[1]);
        (*ri)[2] +=  H * ((*ci)[2] * R2 + fdx * dist[2]);
      }
    }
  }

  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    // modifications needed here
    point_type dist = static_cast<point_type>(source) - center;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);

    auto f0 = charge[0], f1 = charge[1], f2 = charge[2];
    auto fdotx = f0*source[0] + f1*source[1] + f2*source[2];

    for (int n=0; n!=P; ++n) {
      for (int m=0; m<=n; ++m) {
        const int nm  = n * (n + 1) + m;
        const int nms = n * (n + 1) / 2 + m;

        M[0][nms] += f0 * Ynm[nm];
        M[1][nms] += f1 * Ynm[nm];
        M[2][nms] += f2 * Ynm[nm];
        M[3][nms] += fdotx * Ynm[nm];
      }
    }
  }

  void M2M(const multipole_type& Msource,
                 multipole_type& Mtarget, const point_type& translation) const {
    for (unsigned i=0; i<4; i++) {
      LaplaceSpherical::M2M(Msource[i], Mtarget[i], translation);
    }
  }

  /** Kernel M2P operation
   * r += Op(M, t) where M is the multipole and r is the result
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] target The target to evaluate the multipole at
   * @param[in,out] result The target's corresponding result to accumulate into
   * @pre M includes the influence of all sources within its box
   */
  void M2P(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    point_type dist = target - center;
    point_type gradient[4]; //   = {0.,0.,0.,0.};
    gradient[0] = 0.;
    gradient[1] = 0.;
    gradient[2] = 0.;
    gradient[3] = 0.;

    point_type cartesian(0);
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi,Ynm,YnmTheta);

    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      result[0] += std::real(M[0][nms] * Ynm[nm]);
      result[1] += std::real(M[1][nms] * Ynm[nm]);
      result[2] += std::real(M[2][nms] * Ynm[nm]);

      real factor = 1. / r * (n+1);

      gradient[0][0] -= std::real(M[0][nms] * Ynm[nm]) * factor;
      gradient[0][1] += std::real(M[0][nms] * YnmTheta[nm]);

      gradient[1][0] -= std::real(M[1][nms] * Ynm[nm]) * factor;
      gradient[1][1] += std::real(M[1][nms] * YnmTheta[nm]);

      gradient[2][0] -= std::real(M[2][nms] * Ynm[nm]) * factor;
      gradient[2][1] += std::real(M[2][nms] * YnmTheta[nm]);

      gradient[3][0] -= std::real(M[3][nms] * Ynm[nm]) * factor;
      gradient[3][1] += std::real(M[3][nms] * YnmTheta[nm]);

      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        result[0] += 2 * std::real(M[0][nms] * Ynm[nm]);
        result[1] += 2 * std::real(M[1][nms] * Ynm[nm]);
        result[2] += 2 * std::real(M[2][nms] * Ynm[nm]);

        gradient[0][0] -= 2 * std::real(M[0][nms] *Ynm[nm]) * factor;
        gradient[0][1] += 2 * std::real(M[0][nms] *YnmTheta[nm]);
        gradient[0][2] += 2 * std::real(M[0][nms] *Ynm[nm] * CI) * m;

        gradient[1][0] -= 2 * std::real(M[1][nms] *Ynm[nm]) * factor;
        gradient[1][1] += 2 * std::real(M[1][nms] *YnmTheta[nm]);
        gradient[1][2] += 2 * std::real(M[1][nms] *Ynm[nm] * CI) * m;

        gradient[2][0] -= 2 * std::real(M[2][nms] *Ynm[nm]) * factor;
        gradient[2][1] += 2 * std::real(M[2][nms] *YnmTheta[nm]);
        gradient[2][2] += 2 * std::real(M[2][nms] *Ynm[nm] * CI) * m;

        gradient[3][0] -= 2 * std::real(M[3][nms] *Ynm[nm]) * factor;
        gradient[3][1] += 2 * std::real(M[3][nms] *YnmTheta[nm]);
        gradient[3][2] += 2 * std::real(M[3][nms] *Ynm[nm] * CI) * m;
      }
    }
    sph2cart(r,theta,phi,gradient[0],cartesian);
    cartesian *= -target[0];
    gradient[0] = cartesian;

    sph2cart(r,theta,phi,gradient[1],cartesian);
    cartesian *= -target[1];
    gradient[1] = cartesian;

    sph2cart(r,theta,phi,gradient[2],cartesian);
    cartesian *= -target[2];
    gradient[2] = cartesian;

    sph2cart(r,theta,phi,gradient[3],cartesian);
    gradient[3] = cartesian;

    result[0] += (gradient[0][0]+gradient[1][0]+gradient[2][0]+gradient[3][0]);
    result[1] += (gradient[0][1]+gradient[1][1]+gradient[2][1]+gradient[3][1]);
    result[2] += (gradient[0][2]+gradient[1][2]+gradient[2][2]+gradient[3][2]);
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
      LaplaceSpherical::L2L(Lsource[i], Ltarget[i], translation);
    }
  }

 /** Kernel L2P operation
   * r += Op(L, t) where L is the local expansion and r is the result
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] target The target of this L2P operation
   * @param[in] result The result to accumulate into
   * @pre L includes the influence of all sources outside its box
   */
  void L2P(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    point_type dist = target - center;
    point_type gradient[4]; //   = {0.,0.,0.,0.};
    gradient[0] = 0.;
    gradient[1] = 0.;
    gradient[2] = 0.;
    gradient[3] = 0.;
    point_type cartesian(0);

    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi,Ynm,YnmTheta);

    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      result[0] += std::real(L[0][nms] * Ynm[nm]);
      result[1] += std::real(L[1][nms] * Ynm[nm]);
      result[2] += std::real(L[2][nms] * Ynm[nm]);

      real factor = 1. / r * n;
      gradient[0][0] += std::real(L[0][nms] * Ynm[nm]) * factor;
      gradient[0][1] += std::real(L[0][nms] * YnmTheta[nm]);

      gradient[1][0] += std::real(L[1][nms] * Ynm[nm]) * factor;
      gradient[1][1] += std::real(L[1][nms] * YnmTheta[nm]);

      gradient[2][0] += std::real(L[2][nms] * Ynm[nm]) * factor;
      gradient[2][1] += std::real(L[2][nms] * YnmTheta[nm]);

      gradient[3][0] += std::real(L[3][nms] * Ynm[nm]) * factor;
      gradient[3][1] += std::real(L[3][nms] * YnmTheta[nm]);

      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        result[0] += 2 * std::real(L[0][nms] * Ynm[nm]);
        result[1] += 2 * std::real(L[1][nms] * Ynm[nm]);
        result[2] += 2 * std::real(L[2][nms] * Ynm[nm]);

        gradient[0][0] += 2 * std::real(L[0][nms] * Ynm[nm]) * factor;
        gradient[0][1] += 2 * std::real(L[0][nms] * YnmTheta[nm]);
        gradient[0][2] += 2 * std::real(L[0][nms] * Ynm[nm] * CI) * m;

        gradient[1][0] += 2 * std::real(L[1][nms] * Ynm[nm]) * factor;
        gradient[1][1] += 2 * std::real(L[1][nms] * YnmTheta[nm]);
        gradient[1][2] += 2 * std::real(L[1][nms] * Ynm[nm] * CI) * m;

        gradient[2][0] += 2 * std::real(L[2][nms] * Ynm[nm]) * factor;
        gradient[2][1] += 2 * std::real(L[2][nms] * YnmTheta[nm]);
        gradient[2][2] += 2 * std::real(L[2][nms] * Ynm[nm] * CI) * m;

        gradient[3][0] += 2 * std::real(L[3][nms] * Ynm[nm]) * factor;
        gradient[3][1] += 2 * std::real(L[3][nms] * YnmTheta[nm]);
        gradient[3][2] += 2 * std::real(L[3][nms] * Ynm[nm] * CI) * m;
      }
    }
    sph2cart(r,theta,phi,gradient[0],cartesian);
    cartesian *= -target[0];
    gradient[0] = cartesian;

    sph2cart(r,theta,phi,gradient[1],cartesian);
    cartesian *= -target[1];
    gradient[1] = cartesian;

    sph2cart(r,theta,phi,gradient[2],cartesian);
    cartesian *= -target[2];
    gradient[2] = cartesian;

    sph2cart(r,theta,phi,gradient[3],cartesian);
    gradient[3] = cartesian;

    result[0] += (gradient[0][0]+gradient[1][0]+gradient[2][0]+gradient[3][0]);
    result[1] += (gradient[0][1]+gradient[1][1]+gradient[2][1]+gradient[3][1]);
    result[2] += (gradient[0][2]+gradient[1][2]+gradient[2][2]+gradient[3][2]);
  }

}; // end class StokesSpherical
