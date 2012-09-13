#pragma once
/** @file HelmholtzKernel.hpp
 * @brief Helmholtz Kernel class to provide Helmholtz FMM operations.
 */

#include "Point.hpp"

#include <vector>

/** @class HelmholtzKernel
 * @brief The Helmholtz kernel is defined as K(pi,pj) = exp(ik|pi-pj|)/(|pi-pj|).
 * This class implements
 */
class HelmHoltzKernel
{
  /* The wavenumber of the Helmholtz kernel */
  double kappa;

  static constexpr int dimension = 3;
  typedef Point point_type;
  typedef std::complex<double> charge_type;
  typedef std::complex<double> range_type;

  typedef S2Function multipole_type;
  typedef S2Function local_type;

  std::string name() {
    return "Helmholtz Kernel";
  }

  multipole_type init_multipole(double box_size) const {
    return multipole_type();
  }
  local_type init_local(double box_size) const {
    return local_type();
  }

  static constexpr bool is_translation_invariant = true;
  /** Rotation-invariant Kernel evalutation flag */
  static constexpr bool is_rotation_invariant = true;
  /** Symmetric Kernel evaulation flag */
  static constexpr bool is_symmetric = true;
  /** Anti-symmetric Kernel evaluation flag */
  static constexpr bool is_antisymmetric = false;

  /* Only the trivial Kernel is both symmetric and anti-symmetric */
  static_assert(!(is_antisymmetric && is_symmetric),
                "Kernel both symmetric and antisymmetric");
  /* If rotationally-invariant than must also be symmetric */
  static_assert(!is_rotation_invariant || is_symmetric,
                "Kernel rotationally-invariant, but not symmetric");

  /** Translation and rotation-invariant P2P Kernel evaluation
   *
   * @param[in] rij Positive scalar to evaluate the kernel, @a pij = |pi-pj|
   * @result The value of the Kernel: K(pij) = K(|pi-pj|)
   */
  inline range_type operator()(const double pij) const {
    return 1 / pij;
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
