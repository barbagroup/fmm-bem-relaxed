#pragma once
/** @file Kernel.hpp
 * @brief Example Kernel class to be used for Tree/FMM codes.
 */

#include "Point.hpp"

#include <vector>

/** @struct Kernel
 * @brief Example Kernel that can be followed to develop Kernel classes for
 * Tree/FMM code.
 *
 * This class acts as a starting point for defining Kernels that may be
 * used with the treecode and fast multipole method implementations provided in
 * TODO.txt.
 *
 * Not all methods must be implemented, but most constexr and typedefs should
 * be kept for compatibility.
 */
struct Kernel
{
  // Only dimension = 3 support. TODO
  static constexpr int dimension = 3;
  // Spacial point type
  typedef Point point_type;
  // The Kernel return type
  typedef double kernel_value_type;
  // The charge type
  typedef double charge_type;
  // The range type, kernel_value_type * charge_type
  typedef double range_type;

  // The multipole expansion type
  typedef std::vector<double> multipole_type;
  // The local expansion type
  typedef std::vector<double> local_type;

  /** A name identifier for this Kernel */
  std::string name() const {
    return "Kernel Name";
  }

  // TODO
  multipole_type init_multipole(double box_size) const {
    return multipole_type();
  }
  local_type init_local(double box_size) const {
    return local_type();
  }

  /** Kernel evaluation, K(t,s)
   *
   * @param[in] t The target point
   * @param[in] s The source point
   * @result The value of the Kernel, K(t,s)
   */
  inline kernel_value_type operator()(const point_type& t,
                                      const point_type& s) const {
    return 1;
  }


  /** P2P symmetric source-target Kernel accumulation
   * ri += K(pi,pj) * cj
   * rj += K(pj,pi) * ci
   *
   * @param[in]     pi,pj The interacting points
   * @param[in]     ci,cj The charges of pi and pj respectively
   * @param[in,out] ri,rj The results to accumulate into respectively
   */
  inline void P2P(const point_type& pi, const charge_type& ci,
                  range_type& ri,
                  const point_type& pj, const charge_type& cj,
                  range_type& rj) const {
    range_type Kij = operator()(pi,pj);
    ri += Kij * cj;
    rj += Kij * ci;
  }

  /** P2P source-target Kernel accumulation
   * r += K(t,s) * c
   *
   * @param[in]     t The target point
   * @param[in]     s The source point
   * @param[in]     c The source charge
   * @param[in,out] r The result to accumulate into
   */
  inline void P2P(const point_type& t, const point_type& s,
                  const charge_type& c,
                  range_type& r) const {
    r += operator()(t,s) * c;
  }

  /** m += P2M(p,c,r) */
  inline void P2M(const point_type& p, const charge_type& c,
                  const point_type& box_center,
                  multipole_type& m) const {
    return;
  }

  inline void M2M(const multipole_type& m1,
                  const point_type& r,
                  multipole_type& m2) const {
    return;
  }

  inline void M2L(const multipole_type& m,
                  const point_type& r0,
                  local_type& l) const {
    return;
  }

  inline void L2L(const local_type& l1,
                  const point_type& r,
                  local_type& l2) const {
    return;
  }

  inline void L2P(const local_type& l,
                  const point_type& box_center,
                  const point_type& p,
                  range_type& r) const {
    return;
  }

  inline void M2P(const multipole_type& m,
                  const point_type& box_center,
                  const point_type& p,
                  range_type& r) const {
    return;
  }
};
