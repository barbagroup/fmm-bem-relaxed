#pragma once

/** @file StokesSpherical.hpp
 * @brief Implements Stokes kernel with repeated Laplace calls
 */

#include "LaplaceSpherical.hpp"
#include "Mat3.hpp"
#include <iostream>

class StokesSpherical : public LaplaceSpherical
{
 public:
  //! point type
  typedef LaplaceSpherical::point_type point_type;
  //! source type
  typedef LaplaceSpherical::source_type source_type;
  //! target type
  typedef LaplaceSpherical::target_type target_type;
#ifdef STRESSLET
  //! charge type - { g1, g2, g3, n1, n2, n3 }
  typedef Vec<6, LaplaceSpherical::charge_type> charge_type;
#else
  //! charge type - { f1, f2, f3 }
  typedef Vec<3, LaplaceSpherical::charge_type> charge_type;
#endif
  //! return of kernel evaluation
  typedef Mat3<real> kernel_value_type;
  //! product of charge and KVT
  typedef Vec<3, real> result_type;

  //! multipole type
  typedef std::vector<LaplaceSpherical::multipole_type> multipole_type;
  //! local type
  typedef std::vector<LaplaceSpherical::local_type> local_type;

  //! default (delegating) constructor
  StokesSpherical() : StokesSpherical(5) {};
  //! Constructor
  StokesSpherical(int p) : LaplaceSpherical(p) {
#ifdef STRESSLET
    std::cout << "Stresslet calculation" << std::endl;
#endif
  };

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


#ifndef STRESSLET
  /** Kernel evaluation
   * K(t,s)
   * this is kind of useless -- I need to define a P2P anyway, as charge * KVT doesn't work for this app...
   */
  kernel_value_type operator()(const target_type& t, const source_type& s) const {
    // do shit here
    point_type dist = s - t;
    real r2 = normSq(dist);
    real invR2 = 1.0 / r2;
    if (r2 < 1e-8) invR2 = 0;
    real invR3 = invR2 * std::sqrt(invR2);

    kernel_value_type r(0.);
    real dx = dist[0], dy = dist[1], dz = dist[2];

    r(0,0) = invR3*(r2+dx*dx); r(0,1) = invR3*dx*dy; r(0,2) = invR3*dx*dz;
    r(1,0) = invR3*dx*dy; r(1,1) = invR3*(r2+dy*dy); r(1,2) = invR3*dy*dz;
    r(2,0) = invR3*dx*dz; r(2,1) = invR3*dy*dz; r(2,2) = invR3*(r2+dz*dz);
    // real invR = std::sqrt(invR2);
    return r;
  }
#else
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
        real r2 = normSq(dist);
        real invR = 1. / r2;
        if (r2 < 1e-8) invR = 0;

        real dxdotn = dist[0]*(*ci)[3] + dist[1]*(*ci)[4] + dist[2]*(*ci)[5];

        real H = std::sqrt(invR) * invR; // 1 / R^3
        H *= dxdotn * invR;  // (dx . n) / r^5

        real dx0 = dist[0], dx1 = dist[1], dx2 = dist[2];
        auto& g = *ci;

        (*ri)[0] += H * (dx0*dx0*g[0] + dx0*dx1*g[1] + dx0*dx2*g[2]);
        (*ri)[1] += H * (dx0*dx1*g[0] + dx1*dx1*g[1] + dx1*dx2*g[2]);
        (*ri)[2] += H * (dx0*dx2*g[0] + dx1*dx2*g[1] + dx2*dx2*g[2]);
      }
    }
  }
#endif

#ifndef STRESSLET
  /**
   * Create expansions for S_ij / F_i (Tornberg & Greengard
   */
  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    // modifications needed here
    point_type dist = static_cast<point_type>(source) - center;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);

    real f0 = charge[0], f1 = charge[1], f2 = charge[2];
    real fdotx = f0*source[0] + f1*source[1] + f2*source[2];

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
#else
  /**
   * Create expansions for D_ij / G_i (Tornberg & Greengard
   */
  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    // modifications needed here
    point_type dist = static_cast<point_type>(source) - center;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);

    real g0 = charge[0], g1 = charge[1], g2 = charge[2];
    real n0 = charge[3], n1 = charge[4], n2 = charge[5];

    for (int n=0; n!=P; ++n) {
      for (int m=0; m<=n; ++m) {
        const int nm  = n * (n + 1) + m;
        const int nms = n * (n + 1) / 2 + m;

        complex brh = (double)n/rho*Ynm[nm]; // d(rho)
        complex bal = YnmTheta[nm];          // d(alpha)
        complex bbe = -complex(0,1.)*(double)m*Ynm[nm]; // d(beta)

        complex bxd = sin(alpha)*cos(beta)*brh + cos(alpha)*cos(beta)/rho*bal - sin(beta)/rho/sin(alpha)*bbe; // dx
        complex byd = sin(alpha)*sin(beta)*brh + cos(alpha)*sin(beta)/rho*bal + cos(beta)/rho/sin(alpha)*bbe; // dy
        complex bzd = cos(alpha)*brh - sin(alpha)/rho*bal; // dz

        // which order should these be in?
        real rdotn = bxd*n0 + byd*n1 + bzd*n2;
        real rdotg = bxd*g0 + byd*g1 + bzd*g2;
        M[0][nms] += (rdotn * g0 + rdotg * n0);
        M[1][nms] += (rdotn * g1 + rdotg * n1);
        M[2][nms] += (rdotn * g2 + rdotg * n2);

        real xdotg = source[0]*g0 + source[1]*g1 + source[2]*g2;
        real ndotx = n0*source[0] + n1*source[1] + n2*source[2];
        M[3][nms] += rdotn * xdotg + rdotg * ndotx;
      }
    }
  }
#endif
  void M2M(const multipole_type& Msource,
                 multipole_type& Mtarget, const point_type& translation) const {
    LaplaceSpherical::M2M(Msource[0], Mtarget[0], translation);
    LaplaceSpherical::M2M(Msource[1], Mtarget[1], translation);
    LaplaceSpherical::M2M(Msource[2], Mtarget[2], translation);
    LaplaceSpherical::M2M(Msource[3], Mtarget[3], translation);
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
    gradient[0] = point_type(0.);
    gradient[1] = point_type(0.);
    gradient[2] = point_type(0.);
    gradient[3] = point_type(0.);

    point_type cartesian(0);
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi,Ynm,YnmTheta);

#ifdef STRESSLET
    double scale = 1./6;
#else
    double scale = 1.;
#endif

    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      result[0] += scale*std::real(M[0][nms] * Ynm[nm]);
      result[1] += scale*std::real(M[1][nms] * Ynm[nm]);
      result[2] += scale*std::real(M[2][nms] * Ynm[nm]);

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
        result[0] += scale * 2 * std::real(M[0][nms] * Ynm[nm]);
        result[1] += scale * 2 * std::real(M[1][nms] * Ynm[nm]);
        result[2] += scale * 2 * std::real(M[2][nms] * Ynm[nm]);

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

    result[0] += scale*(gradient[0][0]+gradient[1][0]+gradient[2][0]+gradient[3][0]);
    result[1] += scale*(gradient[0][1]+gradient[1][1]+gradient[2][1]+gradient[3][1]);
    result[2] += scale*(gradient[0][2]+gradient[1][2]+gradient[2][2]+gradient[3][2]);
  }

  void M2L(const multipole_type& Msource,
                 local_type& Ltarget, const point_type& translation) const {
    LaplaceSpherical::M2L(Msource[0], Ltarget[0], translation);
    LaplaceSpherical::M2L(Msource[1], Ltarget[1], translation);
    LaplaceSpherical::M2L(Msource[2], Ltarget[2], translation);
    LaplaceSpherical::M2L(Msource[3], Ltarget[3], translation);
  }

  void L2L(const local_type& Lsource,
                 local_type& Ltarget, const point_type& translation) const {
    LaplaceSpherical::L2L(Lsource[0], Ltarget[0], translation);
    LaplaceSpherical::L2L(Lsource[1], Ltarget[1], translation);
    LaplaceSpherical::L2L(Lsource[2], Ltarget[2], translation);
    LaplaceSpherical::L2L(Lsource[3], Ltarget[3], translation);
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
    gradient[0] = point_type(0.);
    gradient[1] = point_type(0.);
    gradient[2] = point_type(0.);
    gradient[3] = point_type(0.);
    point_type cartesian(0);

    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi,Ynm,YnmTheta);

#ifdef STRESSLET
    double scale = 1./6;
#else
    double scale = 1.;
#endif

    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      result[0] += scale*std::real(L[0][nms] * Ynm[nm]);
      result[1] += scale*std::real(L[1][nms] * Ynm[nm]);
      result[2] += scale*std::real(L[2][nms] * Ynm[nm]);

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
        result[0] += scale * 2 * std::real(L[0][nms] * Ynm[nm]);
        result[1] += scale * 2 * std::real(L[1][nms] * Ynm[nm]);
        result[2] += scale * 2 * std::real(L[2][nms] * Ynm[nm]);

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

    result[0] += scale*(gradient[0][0]+gradient[1][0]+gradient[2][0]+gradient[3][0]);
    result[1] += scale*(gradient[0][1]+gradient[1][1]+gradient[2][1]+gradient[3][1]);
    result[2] += scale*(gradient[0][2]+gradient[1][2]+gradient[2][2]+gradient[3][2]);
  }

}; // end class StokesSpherical
