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

#include <complex>
#include <vector>
#include <Vec.hpp>

class CartesianLaplaceKernel 
{
 private:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  const int P;
  //! Epsilon
  static constexpr real EPS = 1e-6;
  //! Imaginary unit
  const complex CI = complex(0,1);
  //! ODD-Even helper function
  inline int ODDEVEN(int n) const {
    return (((n & 1) == 1) ? -1 : 1);
  };

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
  typedef real kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef real result_type;

  //! default constructor - use delegating constructor
  CartesianLaplaceKernel() : CartesianLaplaceKernel(3) {};
  //! Constructor
  CartesianLaplaceKernel(int p)
      : P(p) { };

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, double box_size) const {
    (void) box_size;
    M = std::vector<real>(P*(P+1)*(P+2)/6, 0);
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
    return kernel_value_type(invR);
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
    for( ; t_begin!=t_end; ++t_begin, ++r_begin) {    // Loop over target bodies
      result_type R(0);

      ChargeIter c = c_begin;
      for(auto s = s_begin ; s!=s_end; ++s, ++c) {  //  Loop over source bodies
        auto dist = *t_begin - *s;                                   //   Distance vector from source to target
        real R2 = normSq(dist);                                     //   R^2
        real invR2 = 1.0 / R2;                                    //   1 / R^2
        if( R2 < 1e-8 ) invR2 = 0;                                  //   Exclude self interaction
        real invR = (*c) * std::sqrt(invR2);                   //   potential
        R += invR;                                               //   accumulate potential
      }                                                           //  End loop over source bodies

      (*r_begin) += R;                                           //  potential
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
    for ( ; p_begin!=p_end; ++p_begin, ++c_begin) {
      auto dist = center - *p_begin;
      auto scal = *c_begin;
      M[0] += scal;
      M[1] += scal * dist[0];
      M[2] += scal * dist[1];
      M[3] += scal * dist[2];
      M[4] += scal * dist[0] * dist[0] / 2;
      M[5] += scal * dist[1] * dist[1] / 2;
      M[6] += scal * dist[2] * dist[2] / 2;
      M[7] += scal * dist[0] * dist[1];
      M[8] += scal * dist[1] * dist[2];
      M[9] += scal * dist[2] * dist[0];
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
    const auto dist = translation;
    Mtarget[0] += Msource[0];
    Mtarget[1] += Msource[1] +  dist[0] * Msource[0];
    Mtarget[2] += Msource[2] +  dist[1] * Msource[0];
    Mtarget[3] += Msource[3] +  dist[2] * Msource[0];
    Mtarget[4] += Msource[4] +  dist[0] * Msource[1] + dist[0] * dist[0]  * Msource[0] / 2;
    Mtarget[5] += Msource[5] +  dist[1] * Msource[2] + dist[1] * dist[1]  * Msource[0] / 2;
    Mtarget[6] += Msource[6] +  dist[2] * Msource[3] + dist[2] * dist[2]  * Msource[0] / 2;
    Mtarget[7] += Msource[7] + (dist[0] * Msource[2] + dist[1] * Msource[1] + dist[0] * dist[1] * Msource[0]) / 2;
    Mtarget[8] += Msource[8] + (dist[1] * Msource[3] + dist[2] * Msource[2] + dist[1] * dist[2] * Msource[0]) / 2;
    Mtarget[9] += Msource[9] + (dist[2] * Msource[1] + dist[0] * Msource[3] + dist[2] * dist[0] * Msource[0]) / 2;
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
    const auto dist = translation;
    real R = std::sqrt(norm(dist));
    real R3 = R * R * R;
    real R5 = R3 * R * R;
    Ltarget[0] += Msource[0] / R;
    Ltarget[0] += Msource[1] * (-dist[0] / R3);
    Ltarget[0] += Msource[2] * (-dist[1] / R3);
    Ltarget[0] += Msource[3] * (-dist[2] / R3);
    Ltarget[0] += Msource[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    Ltarget[0] += Msource[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    Ltarget[0] += Msource[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
    Ltarget[0] += Msource[7] * (3 * dist[0] * dist[1] / R5);
    Ltarget[0] += Msource[8] * (3 * dist[1] * dist[2] / R5);
    Ltarget[0] += Msource[9] * (3 * dist[2] * dist[0] / R5);
    Ltarget[1] += Msource[0] * (-dist[0] / R3);
    Ltarget[1] += Msource[1] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    Ltarget[1] += Msource[2] * (3 * dist[0] * dist[1] / R5);
    Ltarget[1] += Msource[3] * (3 * dist[0] * dist[2] / R5);
    Ltarget[2] += Msource[0] * (-dist[1] / R3);
    Ltarget[2] += Msource[1] * (3 * dist[1] * dist[0] / R5);
    Ltarget[2] += Msource[2] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    Ltarget[2] += Msource[3] * (3 * dist[1] * dist[2] / R5);
    Ltarget[3] += Msource[0] * (-dist[2] / R3);
    Ltarget[3] += Msource[1] * (3 * dist[2] * dist[0] / R5);
    Ltarget[3] += Msource[2] * (3 * dist[2] * dist[1] / R5);
    Ltarget[3] += Msource[3] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
    Ltarget[4] += Msource[0] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    Ltarget[5] += Msource[0] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    Ltarget[6] += Msource[0] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
    Ltarget[7] += Msource[0] * (3 * dist[0] * dist[1] / R5);
    Ltarget[8] += Msource[0] * (3 * dist[1] * dist[2] / R5);
    Ltarget[9] += Msource[0] * (3 * dist[2] * dist[0] / R5);
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
    for( ; t_begin!=t_end; ++t_begin, ++r_begin) {
      auto dist = *t_begin - Mcenter;
      real R = std::sqrt(norm(dist));
      real R3 = R * R * R;
      real R5 = R3 * R * R;
      *r_begin += M[0] / R;
      *r_begin += M[1] * (-dist[0] / R3);
      *r_begin += M[2] * (-dist[1] / R3);
      *r_begin += M[3] * (-dist[2] / R3);
      *r_begin += M[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
      *r_begin += M[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
      *r_begin += M[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
      *r_begin += M[7] * (3 * dist[0] * dist[1] / R5);
      *r_begin += M[8] * (3 * dist[1] * dist[2] / R5);
      *r_begin += M[9] * (3 * dist[2] * dist[0] / R5);
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
    auto dist = translation;
    for( int i=0; i<10; ++i ) {
      Ltarget[i] += Lsource[i];
    }
    Ltarget[0] += Lsource[1] * dist[0];
    Ltarget[0] += Lsource[2] * dist[1];
    Ltarget[0] += Lsource[3] * dist[2];
    Ltarget[0] += Lsource[4] * dist[0] * dist[0] / 2;
    Ltarget[0] += Lsource[5] * dist[1] * dist[1] / 2;
    Ltarget[0] += Lsource[6] * dist[2] * dist[2] / 2;
    Ltarget[0] += Lsource[7] * dist[0] * dist[1];
    Ltarget[0] += Lsource[8] * dist[1] * dist[2];
    Ltarget[0] += Lsource[9] * dist[2] * dist[0];
    Ltarget[1] += Lsource[4] * dist[0] * dist[0] / 2;
    Ltarget[1] += Lsource[7] * dist[0] * dist[1];
    Ltarget[1] += Lsource[9] * dist[0] * dist[2];
    Ltarget[2] += Lsource[7] * dist[1] * dist[0];
    Ltarget[2] += Lsource[5] * dist[1] * dist[1] / 2;
    Ltarget[2] += Lsource[8] * dist[1] * dist[2];
    Ltarget[3] += Lsource[9] * dist[2] * dist[0];
    Ltarget[3] += Lsource[8] * dist[2] * dist[1];
    Ltarget[3] += Lsource[6] * dist[2] * dist[2] / 2;
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
           ResultIter r_begin) const
  {
    for ( ; t_begin!=t_end; ++t_begin, ++r_begin) {
      auto dist = *t_begin - Lcenter;
      *r_begin += L[0];
      *r_begin += L[1] * dist[0];
      *r_begin += L[2] * dist[1];
      *r_begin += L[3] * dist[2];
      *r_begin += L[4] * dist[0] * dist[0] / 2;
      *r_begin += L[5] * dist[1] * dist[1] / 2;
      *r_begin += L[6] * dist[2] * dist[2] / 2;
      *r_begin += L[7] * dist[0] * dist[1];
      *r_begin += L[8] * dist[1] * dist[2];
      *r_begin += L[9] * dist[2] * dist[0];
    }
  }
};
