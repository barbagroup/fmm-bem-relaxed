#pragma once
/** @file LaplaceCartesian.hpp
 * @brief Implements the Laplace kernel with spherical expansions.
 *
 * K(t,s) = 1 / |s-t|        // Laplace potential
 * K(t,s) = (s-t) / |s-t|^3  // Laplace force
 */


#include <complex>
#include <vector>
#include <Vec.hpp>

class LaplaceSpherical
{
 protected:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  int P;
  //! \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  std::vector<real> prefactor;
  //! \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  std::vector<real> Anm;
  //! M2L translation matrix \f$ C_{jn}^{km} \f$
  std::vector<complex> Cnm;

  //! Epsilon
  static constexpr real EPS = 1e-12;
  //! Imaginary unit
  const complex CI = complex(0,1);
  //! ODD-Even helper function
  inline int ODDEVEN(int n) const {
    return (((n & 1) == 1) ? -1 : 1);
  };

  //! Custom multipole type
  struct multipole {
    std::vector<complex> M;
    real RCRIT;
    real RMAX;

    //! Convenience method
    complex& operator[](const int i) {
      return M[i];
    }
    //! Convenience method
    const complex& operator[](const int i) const {
      return M[i];
    }
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
  typedef Vec<4,real> kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef Vec<4,real> result_type;

  //! Multipole expansion type
  typedef multipole multipole_type;
  //! Local expansion type
  typedef std::vector<complex> local_type;

  //! default constructor - use delegating constructor
  LaplaceSpherical() : LaplaceSpherical(5) {};
  //! Constructor
  LaplaceSpherical(int p)
      : P(p), prefactor(4*P*P), Anm(4*P*P), Cnm(P*P*P*P) {
    precompute();
  }

  /**
   * precompute all values
   * assume memory correctly allocated
   */
  void precompute()
  {
    for( int n=0; n!=2*P; ++n ) {                               // Loop over n in Anm
      for( int m=-n; m<=n; ++m ) {                              //  Loop over m in Anm
        int nm = n*n+n+m;                                       //   Index of Anm
        int nabsm = abs(m);                                     //   |m|
        real fnmm = EPS;                                        //   Initialize (n - m)!
        for( int i=1; i<=n-m; ++i ) fnmm *= i;                  //   (n - m)!
        real fnpm = EPS;                                        //   Initialize (n + m)!
        for( int i=1; i<=n+m; ++i ) fnpm *= i;                  //   (n + m)!
        real fnma = 1.0;                                        //   Initialize (n - |m|)!
        for( int i=1; i<=n-nabsm; ++i ) fnma *= i;              //   (n - |m|)!
        real fnpa = 1.0;                                        //   Initialize (n + |m|)!
        for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;              //   (n + |m|)!
        prefactor[nm] = std::sqrt(fnma/fnpa);                   //   sqrt( (n - |m|)! / (n + |m|)! )
        Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);              //   (-1)^n / sqrt( (n + m)! / (n - m)! )
      }                                                         //  End loop over m in Anm
    }                                                           // End loop over n in Anm

    for( int j=0, jk=0, jknm=0; j!=P; ++j ) {                   // Loop over j in Cjknm
      for( int k=-j; k<=j; ++k, ++jk ){                         //  Loop over k in Cjknm
        for( int n=0, nm=0; n!=P; ++n ) {                       //   Loop over n in Cjknm
          for( int m=-n; m<=n; ++m, ++nm, ++jknm ) {            //    Loop over m in Cjknm
            const int jnkm = (j+n)*(j+n)+j+n+m-k;               //     Index C_{j+n}^{m-k}
            Cnm[jknm] = std::pow(CI,real(abs(k-m)-abs(k)-abs(m)))//     Cjknm
                * real(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]) * EPS;
          }                                                     //    End loop over m in Cjknm
        }                                                       //   End loop over n in Cjknm
      }                                                         //  End loop over in k in Cjknm
    }                                                           // End loop over in j in Cjknm
  }

  void set_p(int p)
  {
    P = p;
    prefactor.resize(4*P*P);
    Anm.resize(4*P*P);
    Cnm.resize(P*P*P*P);

    // recompute all precomputed values for new P
    precompute();
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M,
                      const point_type& extents, unsigned level) const {
    (void) level;
    M.M = std::vector<complex>(P*(P+1)/2, 0);
    M.RMAX = 0;
    M.RCRIT = extents[0] / 2;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L,
                  const point_type& extents, unsigned level) const {
    (void) extents;  // Quiet warning
    (void) level;
    L = std::vector<complex>(P*(P+1)/2, 0);
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
    point_type dist = s - t;         //   Vector from target to source
    real R2 = normSq(dist);          //   R^2
    real invR2 = 1.0 / R2;           //   1 / R^2
    if (R2 < 1e-8) invR2 = 0;        //   Exclude self interaction
    real invR = std::sqrt(invR2);    //   potential
    dist *= invR2 * invR;            //   force
    return kernel_value_type(invR, dist[0], dist[1], dist[2]);
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
    return kernel_value_type(kst[0], -kst[1], -kst[2], -kst[3]);
  }

  /** Kernel P2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The point source
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    point_type dist = source - center;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      for( int m=0; m<=n; ++m ) {
        const int nm  = n * n + n + m;
        const int nms = n * (n + 1) / 2 + m;
        M[nms] += charge * Ynm[nm];
      }
    }
    M.RMAX = std::max(M.RMAX, norm(dist));
    M.RCRIT = std::min(M.RCRIT, M.RMAX);
  }

  /** Kernel P2M operation
   * M = sum_i Op(p_i) * c_i where M is the multipole and p_i are the points
   *
   * @param[in] p_begin,p_end Iterator pair to the points in this operation
   * @param[in] c_begin Corresponding charge iterator for the points
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   * @pre M is the result of init_multipole
   */
  template <typename SourceIter, typename ChargeIter>
  void P2M(SourceIter p_begin, SourceIter p_end, ChargeIter c_begin,
           const point_type& center, multipole_type& M) const {
    real Rmax = 0;
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    for ( ; p_begin != p_end; ++p_begin, ++c_begin) {
      point_type dist = *p_begin - center;
      real R = norm(dist);
      if( R > Rmax ) Rmax = R;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        for( int m=0; m<=n; ++m ) {
          const int nm  = n * n + n + m;
          const int nms = n * (n + 1) / 2 + m;
          M[nms] += (*c_begin) * Ynm[nm];
        }
      }
    }
    M.RMAX = Rmax;
    M.RCRIT = std::min(M.RCRIT, M.RMAX);
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
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    real Rmax = Mtarget.RMAX;
    real R = norm(translation) + Msource.RCRIT;
    if (R > Rmax) Rmax = R;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,translation);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        const int jk = j * j + j + k;
        const int jks = j * (j + 1) / 2 + k;
        complex M = 0;
        for( int n=0; n<=j; ++n ) {
          for( int m=-n; m<=std::min(k-1,n); ++m ) {
            if( j-n >= k-m ) {
              const int jnkm  = (j - n) * (j - n) + j - n + k - m;
              const int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              const int nm    = n * n + n + m;
              M += Msource[jnkms] * std::pow(CI,real(m-abs(m))) * Ynm[nm]
              * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
          for( int m=k; m<=n; ++m ) {
            if( j-n >= m-k ) {
              const int jnkm  = (j - n) * (j - n) + j - n + k - m;
              const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              const int nm    = n * n + n + m;
              M += std::conj(Msource[jnkms]) * Ynm[nm]
              * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
        Mtarget[jks] += M * EPS;
      }
    }
    Mtarget.RMAX = Rmax;
    Mtarget.RCRIT = std::min(Mtarget.RCRIT, Mtarget.RMAX);
  }

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
    point_type spherical;
    point_type cartesian;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      result[0] += std::real(M[nms] * Ynm[nm]);
      spherical[0] -= std::real(M[nms] * Ynm[nm]) / r * (n+1);
      spherical[1] += std::real(M[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        result[0] += 2 * std::real(M[nms] * Ynm[nm]);
        spherical[0] -= 2 * std::real(M[nms] *Ynm[nm]) / r * (n+1);
        spherical[1] += 2 * std::real(M[nms] *YnmTheta[nm]);
        spherical[2] += 2 * std::real(M[nms] *Ynm[nm] * CI) * m;
      }
    }
    sph2cart(r,theta,phi,spherical,cartesian);
    result[1] += cartesian[0];
    result[2] += cartesian[1];
    result[3] += cartesian[2];
  }

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
    point_type spherical;
    point_type cartesian;
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

 protected:

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta,
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
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom *= rho;                                              //  rho^m
      real rhon = rhom;                                         //  rho^n
      for( int n=m+1; n!=P; ++n ) {                             //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

  //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real rho, real alpha, real beta,
                 complex* Ynm, complex* YnmTheta) const {
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
    for( int m=0; m!=2*P; ++m ) {                               // Loop over m in Ynm
      complex eim = std::exp(CI * real(m * beta));              //  exp(i * m * beta)
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom /= rho;                                              //  rho^(-m-1)
      real rhon = rhom;                                         //  rho^(-n-1)
      for( int n=m+1; n!=2*P; ++n ) {                           //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon /= rho;                                            //   rho^(-n-1)
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

  /** Cartesian to spherical coordinates
   */
  void cart2sph(real& r, real& theta, real& phi,
                point_type dist = point_type()) const {
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
};
