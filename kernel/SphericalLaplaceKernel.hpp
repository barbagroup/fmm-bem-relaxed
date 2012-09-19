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

#pragma once

#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

// #define NO_POINTERS

class SphericalLaplaceKernel
{
 public:
  // TODO MOVE THIS TO PROTECTED
  vect X0;                                                      //!< Center of root cell
  real R0;                                                      //!< Radius of root cell
  C_iter               Ci0;                                     //!< Begin iterator for target cells
  C_iter               Cj0;                                     //!< Begin iterator for source cells



 private:
  real *factorial;                                              //!< Factorial
  real *prefactor;                                              //!< \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  real *Anm;                                                    //!< \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  complex *Cnm;                                                 //!< M2L translation matrix \f$ C_{jn}^{km} \f$

 public:
  //! Constructor
  SphericalLaplaceKernel()
      : factorial(), prefactor(), Anm(), Cnm(),
        X0(0), R0(0), Ci0(), Cj0() {}
  //! Destructor
  ~SphericalLaplaceKernel() {}
  //! Copy constructor
  SphericalLaplaceKernel(const SphericalLaplaceKernel&)
      :factorial(), prefactor(), Anm(), Cnm(),
       X0(0), R0(0), Ci0(), Cj0() {}
  //! Overload assignment
  SphericalLaplaceKernel &operator=(const SphericalLaplaceKernel) {return *this;}

  //! Set center of root cell
  void setX0(vect x0) {X0 = x0;}
  //! Set radius of root cell
  void setR0(real r0) {R0 = r0;}

  //! Get center of root cell
  vect getX0() const {return X0;}
  //! Get radius of root cell
  real getR0() const {return R0;}

  //! Precalculate M2L translation matrix
  void preCalculation() {
    printf("PreCalculation starting\n");
    const complex I(0.,1.);                                     // Imaginary unit
    factorial = new real  [P];                                  // Factorial
    prefactor = new real  [4*P*P];                              // sqrt( (n - |m|)! / (n + |m|)! )
    Anm       = new real  [4*P*P];                              // (-1)^n / sqrt( (n + m)! / (n - m)! )
    Cnm       = new complex [P*P*P*P];                          // M2L translation matrix Cjknm

    factorial[0] = 1;                                           // Initialize factorial
    for( int n=1; n!=P; ++n ) {                                 // Loop to P
      factorial[n] = factorial[n-1] * n;                        //  n!
    }                                                           // End loop to P

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
            Cnm[jknm] = std::pow(I,real(abs(k-m)-abs(k)-abs(m)))//     Cjknm
                * real(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]) * EPS;
          }                                                     //    End loop over m in Cjknm
        }                                                       //   End loop over n in Cjknm
      }                                                         //  End loop over in k in Cjknm
    }                                                           // End loop over in j in Cjknm
    printf("PreCalculation finished\n");
  }

  void initialize() {}

  void P2P(C_iter Ci, C_iter Cj) const {         // Laplace P2P kernel on CPU
    for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {    // Loop over target bodies
      real P0 = 0;                                                //  Initialize potential
      vect F0 = 0;                                                //  Initialize force
      for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
        vect dist = Bi->X - Bj->X -Xperiodic;                     //   Distance vector from source to target
        real R2 = norm(dist) + EPS2;                              //   R^2
        real invR2 = 1.0 / R2;                                    //   1 / R^2
        if( R2 == 0 ) invR2 = 0;                                  //   Exclude self interaction
        real invR = Bj->SRC * std::sqrt(invR2);                   //   potential
        dist *= invR2 * invR;                                     //   force
        P0 += invR;                                               //   accumulate potential
        F0 += dist;                                               //   accumulate force
      }                                                           //  End loop over source bodies

      Bi->TRG[0] += P0;                                           //  potential
      Bi->TRG[1] -= F0[0];                                        //  x component of force
      Bi->TRG[2] -= F0[1];                                        //  y component of force
      Bi->TRG[3] -= F0[2];                                        //  z component of force
      // printf("setting: %lg\n",Bi->TRG[0]);
    }                                                             // End loop over target bodies
  }

  void P2M(C_iter Cj) {
    real Rmax = 0;
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    for( B_iter B=Cj->LEAF; B!=Cj->LEAF+Cj->NCLEAF; ++B ) {
      vect dist = B->X - Cj->X;
      real R = std::sqrt(norm(dist));
      if( R > Rmax ) Rmax = R;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        for( int m=0; m<=n; ++m ) {
          const int nm  = n * n + n + m;
          const int nms = n * (n + 1) / 2 + m;
          Cj->M[nms] += B->SRC * Ynm[nm];
        }
      }
    }
    Cj->RMAX = Rmax;
    Cj->RCRIT = std::min(Cj->R,Rmax);
  }

  void M2M(C_iter Ci) {
    const complex I(0.,1.);
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    real Rmax = Ci->RMAX;
    for( C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; ++Cj ) {
      vect dist = Ci->X - Cj->X;
      real R = std::sqrt(norm(dist)) + Cj->RCRIT;
      if( R > Rmax ) Rmax = R;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
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
                M += Cj->M[jnkms] * std::pow(I,real(m-abs(m))) * Ynm[nm]
                    * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
              }
            }
            for( int m=k; m<=n; ++m ) {
              if( j-n >= m-k ) {
                const int jnkm  = (j - n) * (j - n) + j - n + k - m;
                const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                const int nm    = n * n + n + m;
                M += std::conj(Cj->M[jnkms]) * Ynm[nm]
                    * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
              }
            }
          }
          Ci->M[jks] += M * EPS;
        }
      }
    }
    Ci->RMAX = Rmax;
    Ci->RCRIT = std::min(Ci->R,Rmax);
  }

  void M2L(C_iter Ci, C_iter Cj) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    vect dist = Ci->X - Cj->X - Xperiodic;
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
            L += std::conj(Cj->M[nms]) * Cnm[jknm] * Ynm[jnkm];
          }
          for( int m=0; m<=n; ++m ) {
            const int nm   = n * n + n + m;
            const int nms  = n * (n + 1) / 2 + m;
            const int jknm = jk * P * P + nm;
            const int jnkm = (j + n) * (j + n) + j + n + m - k;
            L += Cj->M[nms] * Cnm[jknm] * Ynm[jnkm];
          }
        }
        Ci->L[jks] += L;
      }
    }
  }

  void M2P(C_iter Ci, C_iter Cj) const {
    const complex I(0.,1.);                                       // Imaginary unit
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
      vect dist = B->X - Cj->X - Xperiodic;
      vect spherical = 0;
      vect cartesian = 0;
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalLocal(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        B->TRG[0] += std::real(Cj->M[nms] * Ynm[nm]);
        spherical[0] -= std::real(Cj->M[nms] * Ynm[nm]) / r * (n+1);
        spherical[1] += std::real(Cj->M[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          B->TRG[0] += 2 * std::real(Cj->M[nms] * Ynm[nm]);
          spherical[0] -= 2 * std::real(Cj->M[nms] *Ynm[nm]) / r * (n+1);
          spherical[1] += 2 * std::real(Cj->M[nms] *YnmTheta[nm]);
          spherical[2] += 2 * std::real(Cj->M[nms] *Ynm[nm] * I) * m;
        }
      }
      sph2cart(r,theta,phi,spherical,cartesian);
      B->TRG[1] += cartesian[0];
      B->TRG[2] += cartesian[1];
      B->TRG[3] += cartesian[2];
    }
  }

  void L2L(C_iter Ci) const {
    const complex I(0.,1.);
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    C_iter Cj = Ci0 + Ci->PARENT;
    vect dist = Ci->X - Cj->X;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
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
            L += std::conj(Cj->L[nms]) * Ynm[jnkm]
                * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
          }
          for( int m=0; m<=n; ++m ) {
            if( n-j >= abs(m-k) ) {
              const int jnkm = (n - j) * (n - j) + n - j + m - k;
              const int nm   = n * n + n + m;
              const int nms  = n * (n + 1) / 2 + m;
              L += Cj->L[nms] * std::pow(I,real(m-k-abs(m-k)))
                  * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
            }
          }
        }
        Ci->L[jks] += L * EPS;
      }
    }
  }

  void L2P(C_iter Ci) const {                                   //!< Evaluate L2P kernel on CPU
    const complex I(0.,1.);                                     // Imaginary unit
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NCLEAF; ++B ) {
      vect dist = B->X - Ci->X;
      vect spherical = 0;
      vect cartesian = 0;
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalMultipole(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        B->TRG[0] += std::real(Ci->L[nms] * Ynm[nm]);
        spherical[0] += std::real(Ci->L[nms] * Ynm[nm]) / r * n;
        spherical[1] += std::real(Ci->L[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          B->TRG[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]);
          spherical[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]) / r * n;
          spherical[1] += 2 * std::real(Ci->L[nms] * YnmTheta[nm]);
          spherical[2] += 2 * std::real(Ci->L[nms] * Ynm[nm] * I) * m;
        }
      }
      sph2cart(r,theta,phi,spherical,cartesian);
      B->TRG[1] += cartesian[0];
      B->TRG[2] += cartesian[1];
      B->TRG[3] += cartesian[2];
    }
  }

  void finalize() {}                                            //!< Finalize kernels

  static int multipole_size(const int level=0) {
    (void) level;  // Quiet warning
    return P*(P+1)/2;
  }
  static int local_size(const int level=0) {
    (void) level;  // Quiet warning
    return P*(P+1)/2;
  }

  //! Free temporary allocations
  void postCalculation() {
    delete[] factorial;                                         // Free factorial
    printf("deleting prefactor\n");
    delete[] prefactor;                                         // Free sqrt( (n - |m|)! / (n + |m|)! )
    printf("done\n");
    delete[] Anm;                                               // Free (-1)^n / sqrt( (n + m)! / (n - m)! )
    delete[] Cnm;                                               // Free M2L translation matrix Cjknm
  }


 private:

  real getBmax(vect const&X, C_iter C) const {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

  void setCenter(C_iter C) const {
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += std::abs(B->SRC);
      X += B->X * std::abs(B->SRC);
    }
    for( C_iter c=Cj0+C->CHILD; c!=Cj0+C->CHILD+C->NCHILD; ++c ) {
      m += std::abs(c->M[0]);
      X += c->X * std::abs(c->M[0]);
    }
    X /= m;
    C->R = getBmax(X,C);
    C->X = X;
  }

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
    const complex I(0.,1.);                                     // Imaginary unit
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1;                                              // Initialize rho^m
    for( int m=0; m!=P; ++m ) {                                 // Loop over m in Ynm
      complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
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
    const complex I(0.,1.);                                     // Imaginary unit
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
    for( int m=0; m!=2*P; ++m ) {                               // Loop over m in Ynm
      complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
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
  void cart2sph(real& r, real& theta, real& phi, vect dist) const {
    r = sqrt(norm(dist))+EPS;                                   // r = sqrt(x^2 + y^2 + z^2) + eps
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
