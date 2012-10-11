#pragma once
/*
  Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/


class SphericalLaplaceKernel
{
 private:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  const int P;
  //! \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  std::vector<real> prefactor;
  //! \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  std::vector<real> Anm;
  //! M2L translation matrix \f$ C_{jn}^{km} \f$
  std::vector<complex> Cnm;

  //! Epsilon
  static constexpr real EPS = 1e-6;
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
  //! Multipole expansion type
  typedef multipole multipole_type;
  //! Local expansion type
  typedef std::vector<complex> local_type;

  //! The dimension of the Kernel
  static constexpr unsigned dimension = 3;
  //! Point type
  // typedef vec<dimension,real> point_type;
  typedef Vec<dimension,real> point_type;
  //! Charge type
  typedef real charge_type;
  //! The return type of a kernel evaluation
  typedef Vec<4,real> kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef Vec<4,real> result_type;

  //! Constructor
  SphericalLaplaceKernel(int p)
      : P(p), prefactor(4*P*P), Anm(4*P*P), Cnm(P*P*P*P) {
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

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, double box_size) const {
    M.M = std::vector<complex>(P*(P+1)/2, 0);
    M.RMAX = 0;
    M.RCRIT = box_size / 2;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, double box_size) const {
    (void) box_size;  // Quiet warning
    L = std::vector<complex>(P*(P+1)/2, 0);
  }

  /** Kernel evaluation
   * K(s,t)
   *
   * @param[in] s,t The source and target points to evaluate the kernel
   */
  kernel_value_type operator()(const point_type& s,
                               const point_type& t) {
    point_type dist = t - s;         //   Distance vector from source to target
    real R2 = normSq(dist);          //   R^2
    real invR2 = 1.0 / R2;           //   1 / R^2
    if( R2 < 1e-8 ) invR2 = 0;       //   Exclude self interaction
    real invR = std::sqrt(invR2);    //   potential
    dist *= invR2 * invR;            //   force
    return kernel_value_type(invR, -dist[0], -dist[1], -dist[2]);
  }

  /** Kernel vectorized non-symmetric P2P operation
   * r_j += sum_i K(s_i,t_j) * c_i
   *
   * @param[in] s_begin,s_end Iterator pair to the source points
   * @param[in] c_begin Iterator to the source charges
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   */
  template <typename PointIter, typename ChargeIter, typename ResultIter>
  void P2P(PointIter s_begin, PointIter s_end, ChargeIter c_begin,
           PointIter t_begin, PointIter t_end, ResultIter r_begin) const
  {
    for( ; t_begin!=t_end; ++t_begin, ++r_begin) {    // Loop over target bodies
      result_type R(0);
      // for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
      ChargeIter c = c_begin;
      for(auto s = s_begin ; s!=s_end; ++s, ++c) {  //  Loop over source bodies
        auto dist = *t_begin - *s;                                   //   Distance vector from source to target
        real R2 = normSq(dist);                                     //   R^2
        real invR2 = 1.0 / R2;                                    //   1 / R^2
        if( R2 < 1e-8 ) invR2 = 0;                                  //   Exclude self interaction
        real invR = (*c) * std::sqrt(invR2);                   //   potential
        dist *= invR2 * invR;                                     //   force
        R[0] += invR;                                               //   accumulate potential
        R[1] += dist[0];                                               //   accumulate force
        R[2] += dist[1];                                               //   accumulate force
        R[3] += dist[2];                                               //   accumulate force
      }                                                           //  End loop over source bodies

      (*r_begin)[0] += R[0];                                           //  potential
      (*r_begin)[1] -= R[1];                                        //  x component of force
      (*r_begin)[2] -= R[2];                                        //  y component of force
      (*r_begin)[3] -= R[3];                                        //  z component of force
      // printf("setting: %lg\n",Bi->TRG[0]);
    }                                                             // End loop over target bodies
  }

  /** Kernel vectorized symmetric P2P operation
   * ...
   */
  template <typename PointIter, typename ChargeIter, typename ResultIter>
  void P2P(PointIter p1_begin, PointIter p1_end, ChargeIter c1_begin,
           PointIter p2_begin, PointIter p2_end, ChargeIter c2_begin,
           ResultIter r1_begin, ResultIter r2_begin) const {
    // TODO...
    (void) p1_begin;
    (void) p1_end;
    (void) c1_begin;
    (void) p2_begin;
    (void) p2_end;
    (void) c2_begin;
    (void) r1_begin;
    (void) r2_begin;
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
  template <typename PointIter, typename ChargeIter>
  void P2M(PointIter p_begin, PointIter p_end, ChargeIter c_begin,
           const point_type& center, multipole_type& M) const {
    real Rmax = 0;
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    // for( B_iter B=C.LEAF; B!=C.LEAF+C.NCLEAF; ++B ) {
    for ( ; p_begin != p_end; ++p_begin, ++c_begin) {
      auto dist = *p_begin - center;
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
   * r_i += Op(M)
   *
   * @param[in] M The multpole expansion
   * @param[in] Mcenter The center of the box with the multipole expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre M includes the influence of all points within its box
   */
  template <typename PointIter, typename ResultIter>
  void M2P(const multipole_type& M, const point_type& Mcenter,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    for( ; t_begin != t_end ; ++t_begin, ++r_begin ) {
      point_type dist = *t_begin - Mcenter;
      point_type spherical(0);
      point_type cartesian(0);
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalLocal(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        (*r_begin)[0] += std::real(M[nms] * Ynm[nm]);
        spherical[0] -= std::real(M[nms] * Ynm[nm]) / r * (n+1);
        spherical[1] += std::real(M[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          (*r_begin)[0] += 2 * std::real(M[nms] * Ynm[nm]);
          spherical[0] -= 2 * std::real(M[nms] *Ynm[nm]) / r * (n+1);
          spherical[1] += 2 * std::real(M[nms] *YnmTheta[nm]);
          spherical[2] += 2 * std::real(M[nms] *Ynm[nm] * CI) * m;
        }
      }
      sph2cart(r,theta,phi,spherical,cartesian);
      (*r_begin)[1] += cartesian[0];
      (*r_begin)[2] += cartesian[1];
      (*r_begin)[3] += cartesian[2];
    }
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Lsource includes the influence of all points outside its box
   */
  void L2L(const local_type& Lsource,
           local_type& Ltarget,
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
            L += std::conj(Lsource[nms]) * Ynm[jnkm]
                * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
          }
          for( int m=0; m<=n; ++m ) {
            if( n-j >= abs(m-k) ) {
              const int jnkm = (n - j) * (n - j) + n - j + m - k;
              const int nm   = n * n + n + m;
              const int nms  = n * (n + 1) / 2 + m;
              L += Lsource[nms] * std::pow(CI,real(m-k-abs(m-k)))
                  * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
            }
          }
        }
        Ltarget[jks] += L * EPS;
      }
    }
  }

  /** Kernel M2P operation
   * r_i += Op(L)
   *
   * @param[in] L The local expansion
   * @param[in] Lcenter The center of the box with the local expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre L includes the influence of all points outside its box
   */
  template <typename PointIter, typename ResultIter>
  void L2P(const local_type& L, const point_type& Lcenter,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];

    for (auto t = t_begin; t != t_end; ++t, ++r_begin) {
      point_type dist = *t - Lcenter;
      point_type spherical(0);
      point_type cartesian(0);
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalMultipole(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        (*r_begin)[0] += std::real(L[nms] * Ynm[nm]);
        spherical[0] += std::real(L[nms] * Ynm[nm]) / r * n;
        spherical[1] += std::real(L[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          (*r_begin)[0] += 2 * std::real(L[nms] * Ynm[nm]);
          spherical[0] += 2 * std::real(L[nms] * Ynm[nm]) / r * n;
          spherical[1] += 2 * std::real(L[nms] * YnmTheta[nm]);
          spherical[2] += 2 * std::real(L[nms] * Ynm[nm] * CI) * m;
        }
      }
      sph2cart(r,theta,phi,spherical,cartesian);
      (*r_begin)[1] += cartesian[0];
      (*r_begin)[2] += cartesian[1];
      (*r_begin)[3] += cartesian[2];
    }
  }

 private:

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
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
  void evalLocal(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
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

  //! Get r,theta,phi from x,y,z
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

  //! Spherical to cartesian coordinates
  template<typename T>
  void sph2cart(real r, real theta, real phi, T spherical, T &cartesian) const {
    cartesian[0] = sin(theta) * cos(phi) * spherical[0]         // x component (not x itself)
        + cos(theta) * cos(phi) / r * spherical[1]
        - sin(phi) / r / sin(theta) * spherical[2];
    cartesian[1] = sin(theta) * sin(phi) * spherical[0]         // y component (not y itself)
        + cos(theta) * sin(phi) / r * spherical[1]
        + cos(phi) / r / sin(theta) * spherical[2];
    cartesian[2] = cos(theta) * spherical[0]                    // z component (not z itself)
        - sin(theta) / r * spherical[1];
  }
};
