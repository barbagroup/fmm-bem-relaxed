#pragma once
/** @file YukawaSpherical.hpp
 * @brief Implements the Laplace kernel with spherical expansions.
 *
 * K(t,s,beta) = e^(beta*|s-t|) / |s-t|        // Laplace potential
 * K(t,s) = (s-t) / |s-t|^3  // Laplace force
 */


#include <complex>
#include <vector>
#include <map>
#include "Vec.hpp"
#include <boost/multi_array.hpp>
#include <tr1/cmath>

// #define USE_FORT_IN

#ifdef USE_FORT_IN
#include "yuk_wrap.h"
#endif

class YukawaSpherical
{
 protected:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  const int P;
  //! Kappa
  const real Kappa;
  //! \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  std::vector<real> prefactor;
  //! \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  std::vector<real> Anm;
  //! M2L translation matrix \f$ C_{jn}^{km} \f$
  std::vector<complex> Cnm;
  //! M2M rotation matrices for +y and -y directions
  typedef boost::multi_array<real,2> Mat2d;
  typedef boost::multi_array<real,3> Mat3d;
  typedef boost::multi_array_types::extent_range range;
  Mat2d sqc;
  Mat3d M2MRotPlus, M2MRotMinus;
  mutable std::map<int,Mat3d> M2MShiftingMatrices;

  //! square root of binomial coefficients
  std::vector<real> C;

  //! Epsilon
  static constexpr real EPS = 1e-6;
  //! Imaginary unit
  const complex CI = complex(0,1);
  //! ODD-Even helper function
  inline int ODDEVEN(int n) const {
    return (((n & 1) == 1) ? -1 : 1);
  };

 public:
  //! The dimension of the Kernel
  static constexpr unsigned dimension = 3;
  //! Point type
  typedef Vec<dimension,real> point_type;
  //! Source type
  typedef point_type source_type;
  //! Target type
  typedef point_type target_type;
  //! Charge type
  typedef real charge_type;
  //! The return type of a kernel evaluation
  typedef real kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef real result_type;

  struct ScaledSeries {
    std::vector<complex> M;
    real scale;
    unsigned level;

    //! Convenience method
    complex& operator[](const int i) {
      return M[i];
    }
    //! Convenience method
    const complex& operator[](const int i) const {
      return M[i];
    }
  };

  //! Multipole expansion type
  typedef ScaledSeries multipole_type;
  //! Local expansion type
  typedef ScaledSeries local_type;

  //! default constructor - use delegating constructor
  YukawaSpherical() : YukawaSpherical(5,0.125) {};
  //! Constructor
  YukawaSpherical(int p, real k)
  : P(p), Kappa(k), prefactor(4*P*P),
    sqc(boost::extents[4*(P+1)][4*(P+1)]),
    M2MRotPlus(boost::extents[P+1][P+1][range(-P,P)]),
    M2MRotMinus(boost::extents[P+1][P+1][range(-P,P)]),
    C(4*P*P) {
    for( int n=0; n!=2*P; ++n ) {
      for( int m=-n; m<=n; ++m ) {
        int nm = n*n+n+m;
        int nabsm = abs(m);
        real fnma = 1.0;
        for( int i=1; i<=n-nabsm; ++i ) fnma *= i;
        real fnpa = 1.0;
        for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;
        prefactor[nm] = ((fnma/fnpa)*(2*n+1));
      }
    }
    // generate sqrt of binomial coefficients
    genBinomialCoefficients(sqc,4*P);

    // generate +/- y rotation matrices (global)
    genRotationMatrices();
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M,
                      const point_type& extent, const unsigned level) const {
    (void) extent;
    M.level = level;
    M.scale = Kappa/pow(2.,level);
    M.M = std::vector<complex>(P*(P+1)/2, 0);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L,
                  const point_type& extent, const unsigned level) const {
    L.level = level;
    L.scale = Kappa*extent[0];
    L.M = std::vector<complex>(P*(P+1)/2, 0);
  }

  /** Kernel evaluation
  * K(t,s)
  *
  * @param[in] t,s The target and source points to evaluate the kernel
  * @result The Laplace potential and force 4-vector on t from s:
  * Potential: 1/|s-t|  Force: (s-t)/|s-t|^3
  */
  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    point_type dist = t - s;
    real r2 = normSq(dist);
    real r  = std::sqrt(r2);
    real invR  = 1.0/r/Kappa;
    if( r < 1e-8 ) invR = 0.;
    real pot = exp(-Kappa*r) * invR * M_PI / 2.;
    return kernel_value_type(pot);
  }

  /** Optional Kernel value source and target transposition
  * K(t,s) -> K(s,t)
  * Often, a kernel has a symmetry in s and t that can be computed faster than
  * by calling the evaluation operator. If this function is implemented, the
  * computation may use it to prevent uneccessary calls to the evaluation
  * operator and accelerate the P2P evaluations and
  *
  * @param[in] kst A kernel value that was returned from operator()(s,t)
  * @returns The value of K(t,s)
  */
  kernel_value_type transpose(const kernel_value_type& kst) const {
    return kernel_value_type(kst);
  }

  /** Kernel P2M operation
  * M += Op(s) * c where M is the multipole and s is the source
  *
  * @param[in] source The point source
  * @param[in] charge The source's corresponding charge
  * @param[in] center The center of the box containing the multipole expansion
  * @param[in,out] M The multipole expansion to accumulate into
  */
  template <typename SourceIter, typename ChargeIter>
  void P2M(SourceIter p_begin, SourceIter p_end, ChargeIter c_begin,
           const point_type& center, multipole_type& M) const {
    real Rmax = 0;
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    double BI[P+1];

    for ( ; p_begin != p_end; ++p_begin, ++c_begin) {
      point_type dist = *p_begin - center;
      real R = norm(dist);
      if( R > Rmax ) Rmax = R;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalLegendre<true>(rho,alpha,-beta,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        for( int m=0; m<=n; ++m ) {
          const int nm  = n * n + n + m;
          const int nms = n * (n + 1) / 2 + m;
          // get the bessel functions
#ifdef USE_FORT_IN
          auto f = Kappa*rho;
          auto Pp1 = P+1;
          double scale = M.scale;
          int calc;
          in_(&scale, &f, &Pp1, BI, &calc);
#else
          In(M.scale, Kappa*rho, P+1, BI);
#endif
          M[nms] += (*c_begin) * Ynm[nm] * BI[n];
        }
      }
    }
  }

  /** Kernel M2M operator
  * M_t += Op(M_s) where M_t is the target and M_s is the source
  *
  * @param[in] source The multipole source at the child level
  * @param[in,out] target The multipole target to accumulate into
  * @param[in] translation The vector from source to target
  * @pre Msource includes the influence of all points within its box
  */
  void M2M(const multipole_type& Msource,
                 multipole_type& Mtarget,
           const point_type& translation) const {
    (void) translation;
    // grab some references
    auto& RD = getRotationMatrix(-translation[2]);
    auto& DC = getShiftingMatrix(Mtarget.level, std::sqrt(normSq(translation)), Msource.scale);
    complex ephi[P+1], Marray[4*P*P];
    const double arg = std::sqrt(2.)/2.;

    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,translation);

    ephi[0] = complex(1.,0);

    // set ephi[1]
    if (translation[0] >= 0) {
      if (translation[1] >= 0) ephi[1] = complex(arg,-arg); // ifl = 3
      else                     ephi[1] = complex(arg,arg);  // ifl = 2
    } else {
      if (translation[1] >= 0) ephi[1] = complex(-arg,-arg);// ifl = 4
      else                     ephi[1] = complex(-arg,arg); // ifl = 1
    }

    // create array of powers of e^{i*m*phi}
    for (int i=1; i<P+1; i++) ephi[i+1] = ephi[i]*ephi[1];

    for (int i=0; i<P+1; i++) {
      // printf("ephi(i): (%.5lg, %.5lg)\n",ephi[i].real(),ephi[i].imag());
    }
    // rotation of phi radians about z-axis in original basis
    multipole_type Mnew;
    init_multipole(Mnew, point_type(0), Mtarget.level);

    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++ ) {
        const int nms = n * (n + 1) / 2 + m;
        Mnew[nms] = conj(ephi[m])*Msource[nms];
      }
    }

    // rotation about y'-axis in rotated system
    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++) {
        const int nm0 = n * (n + 1) / 2;
        const int nms = n * (n + 1) / 2 + m;
        Marray[nms] = Mnew[nm0]*RD[n][0][m];
        for (int mp=1; mp<=n; mp++) {
          const int nmp = n * (n + 1) / 2 + mp;
          Marray[nms] += Mnew[nmp]*RD[n][mp][m]+conj(Mnew[nmp])*RD[n][mp][-m];
        }
      }
    }

    // shift along z-axis
    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++) {
        const int nms = n * (n + 1) / 2 + m;
        Mnew[nms] = 0.;

        for (int mp=m; mp<P; mp++) {
          const int lns = mp * (mp + 1) / 2 + m;
          Mnew[nms] += Marray[lns]*DC[m][n][mp];
        }
      }
    }

    // reverse rotation about y'-axis
    for (int n=0; n<P; n++) {
      // even
      for (int m = 0; m<=n; m+=2) {
        const int nm0 = n * (n + 1) / 2;
        const int nms = n * (n + 1) / 2 + m;
        Marray[nms] = Mnew[nm0]*RD[n][0][m];
        for (int mp = 1; mp <= n; mp += 2) {
          const int nmp = n * (n + 1) / 2 + mp;
          Marray[nms] -= Mnew[nmp]*RD[n][mp][m] + conj(Mnew[nmp])*RD[n][mp][-m];
        }
        for (int mp = 2; mp <= n; mp += 2) {
          const int nmp = n * (n + 1) / 2 + mp;
          Marray[nms] += Mnew[nmp]*RD[n][mp][m] + conj(Mnew[nmp])*RD[n][mp][-m];
        }
      }
      // odd
      for (int m=1; m<=n; m+=2) {
        const int nm0 = n * (n + 1) / 2;
        const int nms = n * (n + 1) / 2 + m;
        Marray[nms] = -Mnew[nm0]*RD[n][0][m];
        for (int mp = 1; mp <= n; mp += 2) {
          const int nmp = n * (n + 1) / 2 + mp;
          Marray[nms] += Mnew[nmp]*RD[n][mp][m] + conj(Mnew[nmp])*RD[n][mp][-m];
        }
        for (int mp=2; mp<=n; mp+=2) {
          const int nmp = n * (n + 1) / 2 + mp;
          Marray[nms] -= Mnew[nmp]*RD[n][mp][m] + conj(Mnew[nmp])*RD[n][mp][-m];
        }
      }
    }

    // rotate back phi radians about z-axis in above system
    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++) {
        const int nms = n * (n + 1) / 2 + m;
        Mnew[nms] = ephi[m]*Marray[nms];
      }
    }

    // append Mnew to Mtarget
    for (unsigned i=0; i<Mnew.M.size(); i++) {
      Mtarget[i] += Mnew[i];
    }
  }

#if 0
  /** Kernel M2L operation
  * L += Op(M)
  *
  * @param[in] Msource The multpole expansion source
  * @param[in,out] Ltarget The local expansion target
  * @param[in] translation The vector from source to target
  * @pre translation obeys the multipole-acceptance criteria
  * @pre Msource includes the influence of all points within its box
  */
  void M2L(const multipole_type& Msource,
           local_type& Ltarget,
      const point_type& translation) const {
      complex Ynm[4*P*P], YnmTheta[4*P*P];

      point_type dist = translation;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalLocal(rho,alpha,beta,Ynm,YnmTheta);
      for( int j=0; j!=P; ++j ) {
        for( int k=0; k<=j; ++k ) {
          const int jk = j * j + j + k;
          const int jks = j * (j + 1) / 2 + k;
          complex L = 0;
          for( int n=0; n!=P; ++n ) {
            for( int m=-n; m<0; ++m ) {
              const int nm   = n * n + n + m;
              const int nms  = n * (n + 1) / 2 - m;
              const int jknm = jk * P * P + nm;
              const int jnkm = (j + n) * (j + n) + j + n + m - k;
              L += std::conj(Msource[nms]) * Cnm[jknm] * Ynm[jnkm];
            }
            for( int m=0; m<=n; ++m ) {
              const int nm   = n * n + n + m;
              const int nms  = n * (n + 1) / 2 + m;
              const int jknm = jk * P * P + nm;
              const int jnkm = (j + n) * (j + n) + j + n + m - k;
              L += Msource[nms] * Cnm[jknm] * Ynm[jnkm];
            }
          }
          Ltarget[jks] += L;
        }
      }
    }
  #endif

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
    double BK[P+1];
    point_type dist = target - center;
    point_type spherical(0);
    point_type cartesian(0);
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLegendre<false>(r,theta,phi,Ynm,YnmTheta);

    // generate k_n(r) terms
    Kn(M.scale,Kappa*r,P+1,BK);

    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      result += std::real(M[nms] * Ynm[nm]) * BK[n];
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        result += 2 * std::real(M[nms] * Ynm[nm]) * BK[n];
      }
    }
  }

  #if 0
  /** Kernel L2L operator
  * L_t += Op(L_s) where L_t is the target and L_s is the source
  *
  * @param[in] source The local source at the parent level
  * @param[in,out] target The local target to accumulate into
  * @param[in] translation The vector from source to target
  * @pre Lsource includes the influence of all points outside its box
  */
  void L2L(const local_type& source,
           local_type& target,
           const point_type& translation) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,translation);
    evalMultipole(rho,alpha,beta,Ynm,YnmTheta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        const int jk = j * j + j + k;
        const int jks = j * (j + 1) / 2 + k;
        complex L = 0;
        for( int n=j; n!=P; ++n ) {
          for( int m=j+k-n; m<0; ++m ) {
            const int jnkm = (n - j) * (n - j) + n - j + m - k;
            const int nm   = n * n + n - m;
            const int nms  = n * (n + 1) / 2 - m;
            L += std::conj(source[nms]) * Ynm[jnkm]
            * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
          }
          for( int m=0; m<=n; ++m ) {
            if( n-j >= abs(m-k) ) {
              const int jnkm = (n - j) * (n - j) + n - j + m - k;
              const int nm   = n * n + n + m;
              const int nms  = n * (n + 1) / 2 + m;
              L += source[nms] * std::pow(CI,real(m-k-abs(m-k)))
              * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
            }
          }
        }
        target[jks] += L * EPS;
      }
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
    point_type spherical(0);
    point_type cartesian(0);
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      result[0] += std::real(L[nms] * Ynm[nm]);
      spherical[0] += std::real(L[nms] * Ynm[nm]) / r * n;
      spherical[1] += std::real(L[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        result[0] += 2 * std::real(L[nms] * Ynm[nm]);
        spherical[0] += 2 * std::real(L[nms] * Ynm[nm]) / r * n;
        spherical[1] += 2 * std::real(L[nms] * YnmTheta[nm]);
        spherical[2] += 2 * std::real(L[nms] * Ynm[nm] * CI) * m;
      }
    }
    sph2cart(r,theta,phi,spherical,cartesian);
    result[1] += cartesian[0];
    result[2] += cartesian[1];
    result[3] += cartesian[2];
  }
#endif

 protected:

  //! Evaluate Legendre Polynomial: Y^m_n(alpha,beta)
  template <bool applyPrefactor>
  void evalLegendre(real rho, real alpha, real beta,
                     complex* Ynm, complex* YnmTheta) const {
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1;                                              // Initialize rho^m
    for( int m=0; m!=P; ++m ) {                                 // Loop over m in Ynm
      complex eim = std::exp(CI * real(m * beta));              //  exp(i * m * beta)
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      if (applyPrefactor) {
        Ynm[npn] = p * prefactor[npn] * eim;               //  rho^m * Ynm for m > 0
      } else {
        Ynm[npn] = p * eim;
      }
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      if (applyPrefactor) {
        YnmTheta[npn] = (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      } else {
        YnmTheta[npn] = (p - (m + 1) * x * p1) / y * eim;
      }
      rhom *= rho;                                              //  rho^m
      real rhon = rhom;                                         //  rho^n
      for( int n=m+1; n!=P; ++n ) {                             //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        if (applyPrefactor) {
          Ynm[npm] = p * prefactor[npm] * eim;             //   rho^n * Ynm
        } else {
          Ynm[npm] = p * eim;
        }
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        if (applyPrefactor) {
          YnmTheta[npm] = ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        } else {
          YnmTheta[npm] = ((n - m + 1) * p - (n + 1) * x * p1) / y * eim;
        }
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

  /** Cartesian to spherical coordinates
  */
  void cart2sph(real& r, real& theta, real& phi, point_type dist=0) const {
    r = norm(dist) + EPS;                                       // r = sqrt(x^2 + y^2 + z^2) + eps
    theta = acos(dist[2] / r);                                  // theta = acos(z / r)
    if( fabs(dist[0]) + fabs(dist[1]) < EPS ) {                 // If |x| < eps & |y| < eps
      phi = 0;                                                  //  phi can be anything so we set it to 0
    } else if( fabs(dist[0]) < EPS ) {                          // If |x| < eps
      phi = dist[1] / fabs(dist[1]) * M_PI * 0.5;               //  phi = sign(y) * pi / 2
    } else if( dist[0] > 0 ) {                                  // If x > 0
      phi = atan(dist[1] / dist[0]);                            //  phi = atan(y / x)
    } else {                                                    // If x < 0
      phi = atan(dist[1] / dist[0]) + M_PI;                     //  phi = atan(y / x) + pi
    }                                                           // End if for x,y cases
  }

  /** Spherical to cartesian coordinates
  */
  template<typename T>
  void sph2cart(real r, real theta, real phi, T spherical,
                T& cartesian) const {
    // x component (not x itself)
    cartesian[0] = sin(theta) * cos(phi) * spherical[0]
                 + cos(theta) * cos(phi) / r * spherical[1]
                 - sin(phi) / r / sin(theta) * spherical[2];
    // y component (not y itself)
    cartesian[1] = sin(theta) * sin(phi) * spherical[0]
                 + cos(theta) * sin(phi) / r * spherical[1]
                 + cos(phi) / r / sin(theta) * spherical[2];
    // z component (not z itself)
    cartesian[2] = cos(theta) * spherical[0]
                 - sin(theta) / r * spherical[1];
  }

  // port of YHROTGEN for only M2M / L2L rotations
  void genRotationMatrices()
  {
    // +y rotation
    double theta = acos(sqrt(3.)/3.);
    yhfstrtn(this->M2MRotPlus, theta);
    // -y rotation
    theta = acos(-sqrt(3.)/3.);
    yhfstrtn(this->M2MRotMinus, theta);
  }

  // get rotation matrix based on +- y
  const Mat3d& getRotationMatrix(double z) const
  {
    if (z < 0.)  return M2MRotPlus;
    else         return M2MRotMinus;
  }

  // port of YHFSTRTN from lapoperators.f:4292
  template <typename Matrix>
  void yhfstrtn(Matrix& D, const real theta)
  {
    // what is this constant?
    double ww = 0.7071067811865476;
    double ctheta = std::cos(theta);
    if (fabs(ctheta) < 1e-8) ctheta = 0.;
    double stheta = std::sin(-theta);
    if (fabs(stheta) < 1e-8) stheta = 0.;

    // more constants
    double hsthta = ww*stheta;
    double cthtap = ww*(1. + ctheta);
    double cthtam = -ww*(1. - ctheta);


    int ij, im, imp;
    // first term is 1.
    D[0][0][0] = 1.;

    for (ij = 1; ij<P; ij++) {

      // m' = 0 case using formula (1)

      for (im = -ij; im <= -1; im++) {
        D[ij][0][im] = -sqc[ij-im][2]*D[ij-1][0][im+1];

        if (im > 1-ij) {
          D[ij][0][im] += sqc[ij+im][2]*D[ij-1][0][im-1];
        }
        D[ij][0][im] *= hsthta;
        if (im > -ij) {
          D[ij][0][im] += D[ij-1][0][im]*ctheta*sqc[ij+im][1]*sqc[ij-im][1];
        }
        D[ij][0][im] /= ij;
      }

      D[ij][0][0] = D[ij-1][0][0]*ctheta;
      if (ij > 1) {
        D[ij][0][0] += hsthta*sqc[ij][2]*( D[ij-1][0][-1] + D[ij-1][0][1] ) / ij;
      }

      for (im = 1; im <= ij; im++) {
        D[ij][0][im] = -sqc[ij+im][2]*D[ij-1][0][im-1];

        if (im < (ij-1)) {
          D[ij][0][im] += sqc[ij-im][2]*D[ij-1][0][im+1];
        }
        D[ij][0][im] *= hsthta;
        if (im < ij) {
          D[ij][0][im] += D[ij-1][0][im]*ctheta*sqc[ij+im][1]*sqc[ij-im][1];
        }
        D[ij][0][im] /= ij;
      }
      // 0 < m' <= j case using formula (2)
      for (imp = 1; imp <= ij; imp++) {
        for (im = -ij; im <= -1; im++) {
          D[ij][imp][im] = D[ij-1][imp-1][im+1]*cthtam*sqc[ij-im][2];

          if (im > (1-ij)) {
            D[ij][imp][im] -= D[ij-1][imp-1][im-1]*cthtap*sqc[ij+im][2];
          }

          if (im > -ij) {
            D[ij][imp][im] += D[ij-1][imp-1][im]*stheta*sqc[ij+im][1]*sqc[ij-im][1];;
          }
          D[ij][imp][im] *= ww/sqc[ij+imp][2];
        }

        D[ij][imp][0] = ij*stheta*D[ij-1][imp-1][0];

        if (ij > 1) {
          D[ij][imp][0] -= sqc[ij][2]*(D[ij-1][imp-1][-1]*cthtap + D[ij-1][imp-1][1]*cthtam);
        }
        D[ij][imp][0] *= ww / sqc[ij+imp][2]; // l: 4406 in lapoperators.f

        for (im = 1; im <= ij; im++) {
          D[ij][imp][im] = D[ij-1][imp-1][im-1]*cthtap*sqc[ij+im][2];

          if (im < (ij-1)) {
            D[ij][imp][im] -= D[ij-1][imp-1][im+1]*cthtam*sqc[ij-im][2];
          }

          if (im < ij) {
            D[ij][imp][im] += D[ij-1][imp-1][im]*stheta*sqc[ij+im][1]*sqc[ij-im][1];
          }
          D[ij][imp][im] *= ww/sqc[ij+imp][2];
        }
      }
    }
    // scale rotation matrix
    double facts[2*P];
    facts[0] = 1.;
    for (int n=1; n<2*P; n++) {
      facts[n] = facts[n-1]*n;
    }

    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++) {
        for (int mp = -n; mp <= n; mp++) {
          auto mpabs = abs(mp);
          D[n][m][mp] *= sqrt(facts[n+m]/facts[n+mpabs]*facts[n-mpabs]/facts[n-m]);
        }
      }
    }
  }

  // compute square root of binomial coefficients
  template <typename Matrix>
  void genBinomialCoefficients(Matrix& C, unsigned nterms)
  {
    // stride of C is P, convert column - row major
    for (unsigned n=0; n<nterms+1; n++) C[n][0] = 1.;

    for (unsigned m=1; m<nterms+1; m++) {
      C[m][m] = 1.;
      for (unsigned n=m+1; n<nterms+1; n++) {
        C[n][m] = C[n-1][m] + C[n-1][m-1];
      }
    }

    // compute the sqrt of these terms
    for (unsigned m=1; m<nterms+1; m++) {
      for (unsigned n=m+1; n<nterms+1; n++) {
        C[n][m] = std::sqrt(C[n][m]);
      }
    }
  }

  // precompute coefficients for shifting of multipole expansion
  // port of YMPSHFTCOEF
  template <typename Matrix>
  Matrix genShiftingMatrix(double scale, double r0) const
  {
    Matrix C0(boost::extents[P+1][P+1][P+1]);
    double fact[2*P+1];
    fact[0] = 1.;
    for (int i=1; i<2*P+1; i++) {
      fact[i] = fact[i-1]*i;
    }

    auto r0k = r0*Kappa;
    double BJ[2*P+2];
#ifdef USE_FORT_IN
    int len = 2*P+2;
    int calc;
    in_(&r0k,&r0k,&len,BJ,&calc);
#else
    In(r0k,r0k,2*P+2,BJ);
#endif

    for (int mnew = 0; mnew <= P; mnew++) {
      for (int lnew = mnew; lnew <= P; lnew++) {
        for (int nn = mnew; nn <= P; nn++) {
          C0[mnew][lnew][nn] = 0.;

          for (int np = mnew; np <= std::min(nn,lnew); np++) {
            C0[mnew][lnew][nn] += pow(scale,nn-lnew)*pow(2.,-lnew-np)*pow(-1.,lnew+nn)*(2*lnew+1)*fact[lnew-mnew]/fact[np+mnew]*fact[nn+mnew] \
                                 *fact[2*np]/fact[np]/fact[np-mnew]/fact[lnew-np]/fact[nn-np]*BJ[lnew+nn-np]*pow(r0k,lnew+nn-2*np);
          }
        }
      }
    }

    return C0;
  }

  // Memoize getting / setting of the shifting matrices
  Mat3d& getShiftingMatrix(int level, double r0, double scale) const
  {
    auto it = M2MShiftingMatrices.find(level);
    if (it != M2MShiftingMatrices.end()) {
      return it->second;
    } else {
      //double diff = 1./pow(2.,level+1) - 1./pow(2.,level);
      //double r0 = std::sqrt(3. * diff*diff);
      //double scale = Kappa/pow(2.,level+1);
      printf("Shifting matrix for level %d, r0 : %.5lg, scale: %.5lg\n",level,r0,scale);
      M2MShiftingMatrices[level].resize(boost::extents[P+1][P+1][P+1]);
      M2MShiftingMatrices[level] = Mat3d(genShiftingMatrix<Mat3d>(scale,r0));
    }

    // if we get here, know shifting matrix has been created
    return M2MShiftingMatrices[level];
  }

  // Modified Spherical Bessel functions
  // direct port of fortran code from Huang et al

  void In(real scale, real x, unsigned nb, real *b) const
  {
    // constants
    real ensig = 1e-4, enmten = 1e-300;

    // if x < 1e-4 then use a 2 term Taylor expansion
    if (x < ensig) {
      real xscale = x / scale;
      real term1 = 1.;
      real term2 = 0.5 * x * x;
      b[0] = term1 * (1. + term2 / 3.);
      for (unsigned i=1; i<=nb; i++) {
        term1 = term1 * xscale / (2*i + 1);

        if (term1 < enmten) term1 = 0.;
        b[i] = term1 * (1. + term2 / (2*i + 3));
      }
    }

    // usually, scal should be less than one, however, in the
    // calculation of h_n, scal will be greater than one,
    // in this case, we should modify the program and replace every
    // small item with zero.
    else if (x > 1e2) {
      for (unsigned i=0; i<=nb; i++) b[i] = 0.;
    } else {
      // otherwise, we will call besselj() and then scale the
      // jn by the scaling factor.
      const real halfpi = M_PI / 2.;
      real cons = std::sqrt(halfpi / x);

      real alpha = 0.5;
      for (unsigned i=0; i<nb+1; i++) {
        b[i] = std::tr1::cyl_bessel_i(alpha+i,x);
      }

      // now scale to i_n*scal**(-n).
      for (unsigned i=0; i<=nb; i++) {
        b[i] *= cons;
        cons /= scale;

        if (fabs(b[i]) < enmten) cons = 0.;
      }
    }
  }

  void Kn(real scale, real x, unsigned nb, real *by) const
  {
    // constants
    const double halfpi = M_PI / 2.;
    const double xmin = 4.46e-308, xinf = 1.79e308, xlarge = 1e8;;
    real ex = x;

    // checks
    if (nb > 0 && x > xmin && x < xlarge) {
      real p = exp(-ex)/ex*halfpi;
      by[0] = p;
      by[1] = p*scale*(1. + 1./ex);

      real u1 = scale / ex;
      real u2 = scale * scale;

      for (unsigned i=2; i<=nb; i++) {
        if (fabs(by[i-1])*u1 > xinf / (2*i-1)) break;

        by[i] = (2*i-1)*u1*by[i-1] + u2*by[i-2];
      }
    } else {
      by[0] = 0.;
    }
  }
};
