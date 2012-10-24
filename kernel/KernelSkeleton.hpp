#pragma once
/*
  Copyright (C) 2012 by TODO

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

/** @class KernelSkeleton
 * @brief Example Kernel that can be followed to develop Kernel classes for
 * Tree/FMM code.
 *
 * This class acts as a starting point for defining Kernels that may be
 * used with the treecode and fast multipole method implementations.
 *
 * Not all methods must be implemented (these are labeled in *Optional*
 * sections), and arbitrary data and helper methods may be added beyond
 * the interface detailed below.
 */
class KernelSkeleton
{
 private:
  //! Precision
  typedef double real;

  // Any other member variables to be used in the Tree operations

 public:
  //! The dimension of the Kernel
  static constexpr unsigned dimension = 3;
  //! Point type -- Can use include/Vec.hpp or anything with op[]
  typedef Vec<dimension,double> point_type;
  //! Charge type
  typedef real charge_type;
  //! The return type of a kernel evaluation
  typedef real kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef real result_type;

  //! Multipole expansion type
  typedef std::vector<real> multipole_type;
  //! Local expansion type
  typedef std::vector<real> local_type;

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, double box_size) const {
    (void) M;
    (void) box_size;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, double box_size) const {
    (void) L;
    (void) box_size;
  }

  /** Kernel evaluation
   * K(t,s)
   *
   * @param[in] t,s The target and source points to evaluate the kernel
   * @return The Kernel evaluation, K(t,s)
   */
  kernel_value_type operator()(const point_type& t,
                               const point_type& s) const {
    (void) t;
    (void) s;
    return kernel_value_type(0);
  }

  /** Kernel P2M operation
   * M = sum_j Op(s_j) * c_j where M is the multipole and s_j are the sources
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
    (void) p_begin;
    (void) p_end;
    (void) c_begin;
    (void) center;
    (void) M;
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Msource includes the influence of all points within its box
   */
  void M2M(const multipole_type& source,
	   multipole_type& target,
           const point_type& translation) const {
    (void) source;
    (void) target;
    (void) translation;
  }

  /** Kernel M2L operation
   * L += Op(M)
   *
   * @param[in] source The multpole expansion source
   * @param[in,out] target The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre Msource includes the influence of all points within its box
   */
  void M2L(const multipole_type& source,
	   local_type& target,
           const point_type& translation) const {
    (void) source;
    (void) target;
    (void) translation;
  }

  /** Kernel M2P operation
   * r_i += Op(M)
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre M includes the influence of all points within its box
   */
  template <typename PointIter, typename ResultIter>
  void M2P(const multipole_type& M, const point_type& center,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    (void) M;
    (void) center;
    (void) t_begin;
    (void) t_end;
    (void) r_begin;
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
    (void) source;
    (void) target;
    (void) translation;
  }

  /** Kernel L2P operation
   * r_i += Op(L)
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre L includes the influence of all points outside its box
   */
  template <typename PointIter, typename ResultIter>
  void L2P(const local_type& L, const point_type& center,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    (void) L;
    (void) center;
    (void) t_begin;
    (void) t_end;
    (void) r_begin;
  }

  /******************/
  /**** Optional ****/
  /******************/
  /* The methods below may be implemented to potentially optimize the P2P operations
   * If these methods are not implemented, the P2P will be delegated to the Direct.hpp
   * methods which use K.operator()(point_type,point_type) for Kernel evaluations.
   */

  /** Kernel vectorized non-symmetric P2P operation
   * r_i += sum_j K(t_i,s_j) * c_j
   *
   * @param[in] s_begin,s_end Iterator pair to the source points
   * @param[in] c_begin Iterator to the source charges
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   */
  template <typename PointIter, typename ChargeIter, typename ResultIter>
  void P2P(PointIter s_begin, PointIter s_end, ChargeIter c_begin,
           PointIter t_begin, PointIter t_end, ResultIter r_begin) const {
    (void) s_begin;
    (void) s_end;
    (void) c_begin;
    (void) t_begin;
    (void) t_end;
    (void) r_begin;
  }

  /** Kernel vectorized symmetric P2P operation
   * r2_i += sum_j K(p2_i, p1_j) * c1_j
   * r1_j += sum_i K(p1_j, p2_i) * c2_i
   *
   * @param[in] p1_begin,p1_end Iterator pair to the source points
   * @param[in] p1_begin Iterator to the source charges
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   */
  template <typename PointIter, typename ChargeIter, typename ResultIter>
  void P2P(PointIter p1_begin, PointIter p1_end, ChargeIter c1_begin,
           PointIter p2_begin, PointIter p2_end, ChargeIter c2_begin,
           ResultIter r1_begin, ResultIter r2_begin) const {
    (void) p1_begin;
    (void) p1_end;
    (void) c1_begin;
    (void) p2_begin;
    (void) p2_end;
    (void) c2_begin;
    (void) r1_begin;
    (void) r2_begin;
  }
};
