#pragma once

#include <complex>
#include <vector>
#include <cassert>
#include <Vec.hpp>

class CartesianYukawaKernel
{
 private:
  typedef double real;
  typedef std::complex<real> complex;
  typedef std::vector<real> real_vec;

  //! Expansion order
  const int P;
  //! Kappa
  const real Kappa;
  //! Epsilon
  static constexpr real EPS = 1e-6;
  //! Imaginary unit
  const complex CI = complex(0,1);
  //! ODD-Even helper function
  inline int ODDEVEN(int n) const {
    return (((n & 1) == 1) ? -1 : 1);
  };
  //! number of multipole terms
  unsigned MTERMS;

  //! I, J, K arrays
  std::vector<unsigned> I, J, K;
  //! indices of multipole terms
  std::vector<unsigned> index;

 private:
  static unsigned setIndex(int P, unsigned i, unsigned j, unsigned k) const {
    unsigned II=0, ii, jj;
    for (unsigned ii = 0; ii < i; ++ii) {
      for (unsigned jj = 1; jj < P+2-ii; ++jj) {
	II += jj;
      }
    }
    for (unsigned jj = P+2-j; jj < P+2; ++jj) {
      II += jj-i;
    }
    return II + k;
  }

 public:
  //! Multipole expansion type
  typedef std::vector<real> multipole_type;
  //! Local expansion type
  typedef std::vector<real> local_type;

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

  //! default constructor - use delegating constructor
  CartesianYukawaKernel() : CartesianYukawaKernel(2,1.) {};
  //! Constructor
  CartesianYukawaKernel(int p, double kappa)
      : P(p), Kappa(kappa), MTERMS((P+1)*(P+2)*(P+3)/6) {
    I = std::vector<unsigned>(MTERMS,0);
    J = std::vector<unsigned>(MTERMS,0);
    K = std::vector<unsigned>(MTERMS,0);
    index = std::vector<unsigned>(MTERMS,0);

    unsigned idx=0;
    for (unsigned ii=0; ii<P+1; ii++) {
      for (unsigned jj=0; jj<P+1-ii; jj++) {
        for (unsigned kk=0; kk<P+1-ii-jj; kk++) {
          index[idx] = setIndex(P,ii,jj,kk);
          I[idx] = ii;
          J[idx] = jj;
          K[idx] = kk;
          idx++;
        }
      }
    }
  };

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, double box_size) const {
    (void) box_size;
    M = std::vector<real>((P+1)*(P+2)*(P+3)/6, 0);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, double box_size) const {
    (void) box_size;  // Quiet warning
    L = std::vector<real>((P+1)*(P+2)*(P+3)/6, 0);
  }

  /** Kernel evaluation
   * K(t,s)
   *
   * @param[in] t,s The target and source points to evaluate the kernel
   * @result The Laplace potential and force 4-vector on t from s:
   * Potential: 1/|s-t|  Force: (s-t)/|s-t|^3
   */
  kernel_value_type operator()(const point_type& t,
                               const point_type& s) {
    point_type dist = s - t;         //   Vector from target to source
    real R2 = normSq(dist);          //   R^2
    real invR2 = 1.0 / R2;           //   1 / R^2
    if (R2 < 1e-8) invR2 = 0;        //   Exclude self interaction
    real invR = std::sqrt(invR2);    //   potential
    dist *= invR2 * invR;            //   force
    return kernel_value_type(invR, dist[0], dist[1], dist[2]);
  }

  /** Kernel vectorized non-symmetric P2P operation
   * r_i += sum_j K(t_i, s_j) * c_j
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
    for( ; t_begin!=t_end; ++t_begin, ++r_begin) {
      result_type R(0);

      ChargeIter c = c_begin;
      for(auto s = s_begin ; s!=s_end; ++s, ++c) {
        auto dist = *t_begin - *s;
        real r2 = normSq(dist);
        real r  = std::sqrt(r2);
        real invR2 = 1.0/r2;
        real invR  = 1.0/r;
        if( r < 1e-8 ) { invR = 0; invR2 = 0; };
        real aux = exp(-Kappa*r)*invR;
        R[0] += (*c)*aux;
        aux *= (Kappa*r+1)*invR2;
        R[1] += (*c)*dist[0]*aux;
        R[2] += (*c)*dist[1]*aux;
        R[3] += (*c)*dist[2]*aux;
      }

      (*r_begin)[0] += R[0];
      (*r_begin)[1] -= R[1];
      (*r_begin)[2] -= R[2];
      (*r_begin)[3] -= R[3];
    }
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
    for ( ; p_begin!=p_end; ++p_begin, ++c_begin) {
      auto dist = center - *p_begin; // dX
      auto scal = *c_begin;

      for (unsigned i=0; i<MTERMS; i++)
      {
        M[i] += scal * pow(dist[0],I[i]) * pow(dist[1],J[i]) * pow(dist[2],K[i]);
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
    const unsigned MTERMS = (P+1)*(P+2)*(P+3)/6;
    const auto dX = translation;

    for (unsigned i=0; i<MTERMS; ++i) {

      for (unsigned ii=0; ii<I[i]+1; ii++) {
        for (unsigned jj=0; jj<J[i]+1; jj++) {
          for (unsigned kk=0; kk<K[i]+1; kk++) {

            unsigned Midx = setIndex(P,ii,jj,kk);
            Mtarget[i] += Msource[Midx]*comb(I[i],I[Midx])*comb(J[i],J[Midx])*comb(K[i],K[Midx])*pow(dX[0],I[i]-I[Midx])*pow(dX[1],J[i]-J[Midx])*pow(dX[2],K[i]-K[Midx]);
          }
        }
      }
    }
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
    (void) Msource;
    (void) Ltarget;
    (void) translation;
    assert(false);
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
           ResultIter r_begin) const
  {

    std::vector<real> a_aux(MTERMS,0), ax_aux(MTERMS,0), ay_aux(MTERMS,0), az_aux(MTERMS,0);

    for( ; t_begin!=t_end; ++t_begin, ++r_begin) {
      auto dX = *t_begin - Mcenter;

      getCoeff(a_aux,ax_aux,ay_aux,az_aux,dX);

      // loop over tarSize
      for (unsigned j=0; j<MTERMS; j++) {
        (*r_begin)[0] +=  a_aux[j]*M[j];
        (*r_begin)[1] += ax_aux[j]*M[j];
        (*r_begin)[2] += ay_aux[j]*M[j];
        (*r_begin)[3] += az_aux[j]*M[j];
      }
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
           const point_type& translation) const
  {
    (void) Lsource;
    (void) Ltarget;
    (void) translation;
    assert(false);
  }

  /** Kernel L2P operation
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
           ResultIter r_begin) const
  {
    (void) L;
    (void) Lcenter;
    (void) t_begin;
    (void) t_end;
    (void) r_begin;
    assert(false);
  }

 private:
  /* No longer needed?
  static constexpr unsigned fact(unsigned N) const {
  // WARNING: Overflows when N > 13
    unsigned f = 1;
    for (unsigned ii = 1; ii < N+1; ++ii)
      f *= ii;
    return f;
  }
  */

  /** Compute "n choose k" combinatorial
   * @returns n! / (k! (n-k)!)
   */
  unsigned comb(unsigned n, unsigned k) const
  {
    // WARNING: Overflows really quickly
    //return fact(n) / (fact(k)*fact(n-k));

    // Better (and faster) way:
    if (k > n) return 0;
    if (2*k > n) k = n-k;
    if (k == 0) return 1;

    unsigned result = n;
    for (unsigned i = 2; i <= k; ++i) {
      result *= (n-i+1);
      result /= i;
    }
    return result;
  }

  void getCoeff(real_vec& a, real_vec& ax, real_vec& ay, real_vec& az, point_type& dX) const
  {
    real_vec b(a.size(),0);
    real dx = dX[0], dy = dX[1], dz = dX[2];

    auto R2 = normSq(dX);
    auto R  = sqrt(R2);

    int i,j,k,I,Im1x,Im2x,Im1y,Im2y,Im1z,Im2z;
    real C,C1,C2,Cb, R2_1;

    R2_1 = 1/R2;

    // First coefficient
    b[0] = exp(-Kappa*R);
    a[0] = b[0]/R;

    // Two indices = 0
    I = setIndex(P,1,0,0);
    b[I]   = -Kappa * (dx*a[0]); // 1,0,0
    b[P+1] = -Kappa * (dy*a[0]); // 0,1,0
    b[1]   = -Kappa * (dz*a[0]); // 0,0,1

    a[I]   = -R2_1*dx*(Kappa*b[0]+a[0]);
    a[P+1] = -R2_1*dy*(Kappa*b[0]+a[0]);
    a[1]   = -R2_1*dz*(Kappa*b[0]+a[0]);

    ax[0]  = a[I];
    ay[0]  = a[P+1];
    az[0]  = a[1];

    for (i=2; i<P+1; i++)
    {
      Cb   = -Kappa/i;
      C    = R2_1/i;
      I    = setIndex(P,i,0,0);
      Im1x = setIndex(P,i-1,0,0);
      Im2x = setIndex(P,i-2,0,0);
      b[I] = Cb * (dx*a[Im1x] + a[Im2x]);
      a[I] = C * ( -Kappa*(dx*b[Im1x] + b[Im2x]) -(2*i-1)*dx*a[Im1x] - (i-1)*a[Im2x] );
      ax[Im1x] = a[I]*i;

      I    = setIndex(P,0,i,0);
      Im1y = I-(P+2-i);
      Im2y = Im1y-(P+2-i+1);
      b[I] = Cb * (dy*a[Im1y] + a[Im2y]);
      a[I] = C * ( -Kappa*(dy*b[Im1y] + b[Im2y]) -(2*i-1)*dy*a[Im1y] - (i-1)*a[Im2y] );
      ay[Im1y] = a[I]*i;

      I   = i;
      Im1z = I-1;
      Im2z = I-2;
      b[I] = Cb * (dz*a[Im1z] + a[Im2z]);
      a[I] = C * ( -Kappa*(dz*b[Im1z] + b[Im2z]) -(2*i-1)*dz*a[Im1z] - (i-1)*a[Im2z] );
      az[Im1z] = a[I]*i;
    }
  // One index = 0, one = 1 other >=1

    Cb   = -Kappa/2;
    C    = R2_1/2.;
    I    = setIndex(P,1,1,0);
    Im1x = P+1;
    Im1y = I-P;
    b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y]);
    a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]) - 3*(dx*a[Im1x]+dy*a[Im1y]) );
    ax[Im1x] = a[I];
    ay[Im1y] = a[I];
    I    = setIndex(P,1,0,1);
    Im1x = 1;
    Im1z = I-1;
    b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z]);
    a[I] = C * ( -Kappa*(dx*b[Im1x]+dz*b[Im1z]) - 3*(dx*a[Im1x]+dz*a[Im1z]) );
    ax[Im1x] = a[I];
    az[Im1z] = a[I];
    I    = setIndex(P,0,1,1);
    Im1y = I-(P+1);
    Im1z = I-1;
    b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z]);
    a[I] = C * ( -Kappa*(dy*b[Im1y]+dz*b[Im1z]) - 3*(dy*a[Im1y]+dz*a[Im1z]) );
    ay[Im1y] = a[I];
    az[Im1z] = a[I];

    for (i=2; i<P; i++)
    {
      Cb   = -Kappa/(i+1);
      C    = R2_1/(1+i);
      C1   = 1+2*i;
      I    = setIndex(P,1,i,0);
      Im1x = setIndex(P,0,i,0);
      Im1y = I-(P+1-i);
      Im2y = Im1y-(P+2-i);
      b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2y]);
      a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
      ax[Im1x] = a[I];
      ay[Im1y] = a[I]*i;

      I    = setIndex(P,1,0,i);
      Im1x = setIndex(P,0,0,i);
      Im1z = I-1;
      Im2z = I-2;
      b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2z]);
      a[I] = C * ( -Kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2z]) - C1*(dx*a[Im1x]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
      ax[Im1x] = a[I];
      az[Im1z] = a[I]*i;

      I    = setIndex(P,0,1,i);
      Im1y = I-(P+1);
      Im1z = I-1;
      Im2z = I-2;
      b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
      a[I] = C * ( -Kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) - C1*(dy*a[Im1y]+dz*a[Im1z]) - (1+i-1)*(a[Im2z]) );
      ay[Im1y] = a[I];
      az[Im1z] = a[I]*i;

      I    = setIndex(P,i,1,0);
      Im1y = I-(P+1-i);
      Im1x = setIndex(P,i-1,1,0);
      Im2x = setIndex(P,i-2,1,0);
      b[I] = Cb * (dy*a[Im1y] + dx*a[Im1x] + a[Im2x]);
      a[I] = C * ( -Kappa*(dy*b[Im1y]+dx*b[Im1x]+b[Im2x]) - C1*(dy*a[Im1y]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
      ax[Im1x] = a[I]*i;
      ay[Im1y] = a[I];

      I    = setIndex(P,i,0,1);
      Im1z = I-1;
      Im1x = setIndex(P,i-1,0,1);
      Im2x = setIndex(P,i-2,0,1);
      b[I] = Cb * (dz*a[Im1z] + dx*a[Im1x] + a[Im2x]);
      a[I] = C * ( -Kappa*(dz*b[Im1z]+dx*b[Im1x]+b[Im2x]) - C1*(dz*a[Im1z]+dx*a[Im1x]) - (1+i-1)*(a[Im2x]) );
      ax[Im1x] = a[I]*i;
      az[Im1z] = a[I];

      I    = setIndex(P,0,i,1);
      Im1z = I-1;
      Im1y = I-(P+2-i);
      Im2y = Im1y-(P+3-i);
      b[I] = Cb * (dz*a[Im1z] + dy*a[Im1y] + a[Im2y]);
      a[I] = C * ( -Kappa*(dz*b[Im1z]+dy*b[Im1y]+b[Im2y]) - C1*(dz*a[Im1z]+dy*a[Im1y]) - (1+i-1)*(a[Im2y]) );
      ay[Im1y] = a[I]*i;
      az[Im1z] = a[I];
    }

    // One index 0, others >=2
    for (i=2; i<P+1; i++)
    {
      for (j=2; j<P+1-i; j++)
      {
        Cb   = -Kappa/(i+j);
        C    = R2_1/(i+j);
        C1   = 2*(i+j)-1;
        I    = setIndex(P,i,j,0);
        Im1x = setIndex(P,i-1,j,0);
        Im2x = setIndex(P,i-2,j,0);
        Im1y = I-(P+2-j-i);
        Im2y = Im1y-(P+3-j-i);
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + a[Im2x] + a[Im2y]);
        a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+b[Im2x]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]) -(i+j-1)*(a[Im2x]+a[Im2y]) );
        ax[Im1x] = a[I]*i;
        ay[Im1y] = a[I]*j;

        I    = setIndex(P,i,0,j);
        Im1x = setIndex(P,i-1,0,j);
        Im2x = setIndex(P,i-2,0,j);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dx*a[Im1x] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
        a[I] = C * ( -Kappa*(dx*b[Im1x]+dz*b[Im1z]+b[Im2x]+b[Im2z]) - C1*(dx*a[Im1x]+dz*a[Im1z]) -(i+j-1)*(a[Im2x]+a[Im2z]) );
        ax[Im1x] = a[I]*i;
        az[Im1z] = a[I]*j;

        I    = setIndex(P,0,i,j);
        Im1y = I-(P+2-i);
        Im2y = Im1y-(P+3-i);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
        a[I] = C * ( -Kappa*(dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) - C1*(dy*a[Im1y]+dz*a[Im1z]) -(i+j-1)*(a[Im2y]+a[Im2z]) );
        ay[Im1y] = a[I]*i;
        az[Im1z] = a[I]*j;
      }
    }

    if (P>2)
    {
      // Two index = 1, other>=1
      C    = R2_1/3;
      Cb   = -Kappa/3;
      I    = setIndex(P,1,1,1);
      Im1x = setIndex(P,0,1,1);
      Im1y = setIndex(P,1,0,1);
      Im1y = I-(P);
      Im1z = I-1;
      b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z]);
      a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]) - 5*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) );
      ax[Im1x] = a[I];
      ay[Im1y] = a[I];
      az[Im1z] = a[I];
      for (i=2; i<P-1; i++)
      {
        Cb   = -Kappa/(2+i);
        C    = R2_1/(i+2);
        C1   = 2*i+3;
        I    = setIndex(P,i,1,1);
        Im1x = setIndex(P,i-1,1,1);
        Im1y = I-(P+1-i);
        Im1z = I-1;
        Im2x = setIndex(P,i-2,1,1);
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x]);
        a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2x]) );
        ax[Im1x] = a[I]*i;
        ay[Im1y] = a[I];
        az[Im1z] = a[I];

        I    = setIndex(P,1,i,1);
        Im1x = setIndex(P,0,i,1);
        Im1y = I-(P+1-i);
        Im2y = Im1y-(P+2-i);
        Im1z = I-1 ;
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y]);
        a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2y]) );
        ax[Im1x] = a[I];
        ay[Im1y] = a[I]*i;
        az[Im1z] = a[I];

        I    = setIndex(P,1,1,i);
        Im1x = setIndex(P,0,1,i);
        Im1y = I-(P);
        Im1z = I-1;
        Im2z = I-2;
        b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2z]);
        a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2z]) - C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - (i+1)*(a[Im2z]) );
        ax[Im1x] = a[I];
        ay[Im1y] = a[I];
        az[Im1z] = a[I]*i;
      }
    }

    // One index = 1, others >=2
    if (P>4)
    {
      for (i=2; i<P-2; i++)
      {
        for (j=2; j<P-i; j++)
        {
          Cb = -Kappa/(1+i+j);
          C  =  R2_1/(1+i+j);
          C1 = -(2.*(i+j)+1);
          C2 = (i+j);
          I    = setIndex(P,1,i,j);
          Im1x = setIndex(P,0,i,j);
          Im1y = I-(P+1-i);
          Im2y = Im1y-(P+2-i);
          Im1z = I-1;
          Im2z = I-2;
          b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2y] + a[Im2z]);
          a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2y]+a[Im2z]) );
          ax[Im1x] = a[I];
          ay[Im1y] = a[I]*i;
          az[Im1z] = a[I]*j;

          I    = setIndex(P,i,1,j);
          Im1x = setIndex(P,i-1,1,j);
          Im1y = I-(P+1-i);
          Im2x = setIndex(P,i-2,1,j);
          Im1z = I-1;
          Im2z = I-2;
          b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2z]);
          a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2z]) );
          ax[Im1x] = a[I]*i;
          ay[Im1y] = a[I];
          az[Im1z] = a[I]*j;

          I    = setIndex(P,i,j,1);
          Im1x = setIndex(P,i-1,j,1);
          Im2x = setIndex(P,i-2,j,1);
          Im1y = I-(P+2-i-j);
          Im2y = Im1y-(P+3-i-j);
          Im1z = I-1;
          b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2y]);
          a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2y]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2y]) );
          ax[Im1x] = a[I]*i;
          ay[Im1y] = a[I]*j;
          az[Im1z] = a[I];
        }
      }
    }

    // All indices >= 2
    if (P>5)
    {
      for (i=2;i<P-3;i++)
      {
        for (j=2;j<P-1-i;j++)
        {
          for (k=2;k<P+1-i-j;k++)
          {
            Cb = -Kappa/(i+j+k);
            C  = R2_1/(i+j+k);
            C1 = -(2.*(i+j+k)-1);
            C2 = i+j+k-1.;
            I    = setIndex(P,i,j,k);
            Im1x = setIndex(P,i-1,j,k);
            Im2x = setIndex(P,i-2,j,k);
            Im1y = I-(P+2-i-j);
            Im2y = Im1y-(P+3-i-j);
            Im1z = I-1;
            Im2z = I-2;
            b[I] = Cb * (dx*a[Im1x] + dy*a[Im1y] + dz*a[Im1z] + a[Im2x] + a[Im2y] + a[Im2z]);
            a[I] = C * ( -Kappa*(dx*b[Im1x]+dy*b[Im1y]+dz*b[Im1z]+b[Im2x]+b[Im2y]+b[Im2z]) + C1*(dx*a[Im1x]+dy*a[Im1y]+dz*a[Im1z]) - C2*(a[Im2x]+a[Im2y]+a[Im2z]) );
            ax[Im1x] = a[I]*i;
            ay[Im1y] = a[I]*j;
            az[Im1z] = a[I]*k;
          }
        }
      }
    }
  }
};
