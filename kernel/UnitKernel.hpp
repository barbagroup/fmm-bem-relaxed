#pragma once
/** @file UnitKernel.hpp
 * @brief Implements the unit kernel defined by
 * K(t,s) = 1  if t != s
 * K(t,s) = 0  if t == s
 */

#include "Vec.hpp"

class UnitKernel
{
 public:
  //! The dimension of the Kernel
  static constexpr unsigned dimension = 3;
  //! Point type
  typedef Vec<dimension,double> point_type;
  //! Source type
  typedef point_type source_type;
  //! Target type
  typedef point_type target_type;
  //! Charge type
  typedef double charge_type;
  //! The return type of a kernel evaluation
  typedef unsigned kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef double result_type;

  //! Multipole expansion type
  typedef double multipole_type;
  //! Local expansion type
  typedef double local_type;

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, const point_type&, unsigned) const {
    M = 0;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, const point_type&, unsigned) const {
    L = 0;
  }

  /** Kernel evaluation
   * K(t,s)
   *
   * @param[in] t,s The target and source points to evaluate the kernel
   */
  kernel_value_type operator()(const point_type& t,
                               const point_type& s) const {
    return t == s ? kernel_value_type(0) : kernel_value_type(1);
  }

  /** Kernel P2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The point source
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void P2M(const point_type&, const charge_type& charge,
           const point_type&, multipole_type& M) const {
    M += charge;
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre source includes the influence of all points within its box
   */
  void M2M(const multipole_type& source,
           multipole_type& target,
           const point_type&) const {
    target += source;
  }

  /** Kernel M2L operation
   * L += Op(M)
   *
   * @param[in] source The multpole expansion source
   * @param[in,out] target The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre source includes the influence of all points within its box
   */
  void M2L(const multipole_type& source,
                 local_type& target,
           const point_type&) const {
    target += source;
  }

  /** Kernel M2P operation
   * r += Op(M) where M is the multipole and r is the result
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] target The target point position
   * @param[in,out] result The target's corresponding result to accumulate into
   * @pre M includes the influence of all points within its box
   */
  void M2P(const multipole_type& M, const point_type&,
           const point_type&, result_type& result) const {
    result += M;
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre source includes the influence of all points outside its box
   */
  void L2L(const local_type& source,
           local_type& target,
           const point_type&) const {
    target += source;
  }

  /** Kernel L2P operation
   * r += Op(L) where L is the local expansion and r is the result
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre L includes the influence of all points outside its box
   */
  void L2P(const local_type& L, const point_type&,
           const point_type&, result_type& result) const {
    result += L;
  }
};
