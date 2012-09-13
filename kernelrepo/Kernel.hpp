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
  static constexpr int dimension = 3;
  typedef Point point_type;
  typedef double charge_type;
  typedef double range_type;

  typedef std::vector<double> multipole_type;
  typedef std::vector<double> local_type;

  std::string name() {
    return "Kernel Name";
  }

  multipole_type init_multipole(double box_size) const {
    return multipole_type();
  }
  local_type init_local(double box_size) const {
    return local_type();
  }

  /** Translation-invariant Kernel flag.
   Translation-variant Kernels are currently not supported.
   Are there any FMMs that do? */
  static constexpr bool is_translation_invariant = true;
  /** Rotation-invariant Kernel evalutation flag */
  static constexpr bool is_rotation_invariant = true;
  /** Symmetric Kernel evaulation flag
   * True when K(r1,r2) = K(r2,r1)
   * This is often the case with potential kernels */
  static constexpr bool is_symmetric = true;
  /** Anti-symmetric Kernel evaluation flag
   * True when K(r1,r2) = -K(r2,r1)
   * This is often the case with force kernels */
  static constexpr bool is_antisymmetric = false;

  /* Only the trivial Kernel is both symmetric and anti-symmetric */
  static_assert(!(is_antisymmetric && is_symmetric),
                "Kernel both symmetric and antisymmetric");
  /* If rotationally-invariant than must also be symmetric */
  static_assert(!is_rotation_invariant || is_symmetric,
                "Kernel rotationally-invariant, but not symmetric");

  /** Translation and rotation-invariant P2P Kernel evaluation
   * Implement this function (and flag) when K(pi,pj) = K(|pi-pj|)
   * This is often the case with potential kernels
   *
   * @param[in] rij Positive scalar to evaluate the kernel, @a pij = |pi-pj|
   * @result The value of the Kernel: K(pij) = K(|pi-pj|)
   */
  inline range_type operator()(const double pij) const {
    return 1;
  }

  /** P2P Kernel evaluation
   *
   * @param[in] pi The position of the ith point
   * @param[in] pj The position of the jth point
   * @result The value of the Kernel: K(pi,pj)
   */
  inline range_type operator()(const point_type& pi,
                               const point_type& pj) const {
    return operator()((pi-pj).magnitude());
  }


  /** Alternative to the above flags and signatures?
   * ::NO::, because want to cut the direct loop in half if symm/antisymm
   */
  inline void P2P(const point_type& pi, const charge_type& ci,
                  const point_type& pj, const charge_type& cj,
                  range_type& ri, range_type& rj) const {
    range_type val = operator()(pi,pj);
    ri += cj * val;
    rj += ci * val;
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
