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


class CountingKernel
{
 private:
  //! Precision
  typedef double real;
  typedef std::complex<real> complex;

 public:
  //! Multipole expansion type
  typedef unsigned multipole_type;
  //! Local expansion type
  typedef unsigned local_type;

  //! The dimension of the Kernel
  static constexpr unsigned dimension = 3;
  //! Point type
  // typedef vec<dimension,real> point_type;
  typedef Vec<dimension,real> point_type;
  //! Charge type
  typedef real charge_type;
  //! The return type of a kernel evaluation
  typedef unsigned kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef unsigned result_type;

  //! Constructor
  CountingKernel() {};

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, double box_size) const {
    M = 0;
    (void) box_size;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, double box_size) const {
    L = 0;
    (void) box_size;
  }

  /** Kernel evaluation
   * K(s,t)
   *
   * @param[in] s,t The source and target points to evaluate the kernel
   */
  kernel_value_type operator()(const point_type& s,
                               const point_type& t) {
    (void) s;
    (void) t;
    return kernel_value_type(1);
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
    (void) c_begin;
    for ( ; t_begin!=t_end; ++t_begin, ++r_begin)
    {
      result_type r(0);
      auto s = s_begin;
      for ( ; s!= s_end; ++s) r++;
      *r_begin += r;
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
           const point_type& center, multipole_type& M) {
    (void) p_begin;
    (void) p_end;
    for ( ; p_begin != p_end; ++p_begin) M++;
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
           const point_type& translation) {
    (void) translation;
    Mtarget += Msource;
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
           const point_type& translation) {
    (void) translation;
    Ltarget += Msource;
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
    (void) Mcenter;
    for ( ; t_begin != t_end; ++t_begin, ++r_begin) *r_begin += M;
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
    (void) translation;
    Ltarget += Lsource;
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
           ResultIter r_begin) const {
    (void) Lcenter;
    for ( ; t_begin!=t_end; ++t_begin, ++r_begin) *r_begin += L;
  }
};
