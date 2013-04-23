#pragma once
/** @file KernelSkeleton
 * @brief An example Kernel implementation that explains the required and
 * optional methods and types that define a Kernel class.
 *
 * This class will be used to evaluate the matrix-vector product
 * r_i = sum_j K(t_i, s_j) c_j
 * where K is the Kernel defined by this class. A kernel is a function:
 * target_type x source_type -> kernel_value_type
 * The types in the matrix-vector product are
 * kernel_value_type x charge_type -> result_type
 *
 * The source_type and target_type must be castable to a point_type. FMMs and
 * Treecodes use spacial trees to partition the sources and targets so they must
 * define a dimension and point_type.
 *
 * Many Kernels will specify source, target, and point types that are the same.
 * Other Kernels must ensure that
 * static_cast<point_type>(target_type)
 * and
 * static_cast<point_type>(source_type)
 * operate as expected.
 *
 * MORE on Tree Operations.
 */


class KernelSkeleton
{
 private:
  //! Precision
  typedef double real;

  // Any other member variables to be used in the Tree operations

 public:
  //! The dimension of the point_type
  static constexpr unsigned dimension = 3;
  //! Point type -- Recommend include/Vec.hpp.
  // TODO Ability to use anything with op[]/size() or begin()/end()?
  typedef Vec<dimension,double> point_type;

  //! Return type of a kernel evaluation
  typedef real kernel_value_type;

  //! Source type
  //! Must be convertible to point_type: static_cast<point_type>(source_type)
  typedef point_type source_type;
  //! Charge type associated with each source
  //! The type of the vector in the FMM/Treecode matvec
  typedef real charge_type;

  //! Target type
  //! Must be convertible to point_type: static_cast<point_type>(target_type)
  typedef point_type target_type;
  //! Result type associated with each target
  //! The product of the kernel_value_type and the charge_type
  typedef real result_type;

  //! Multipole expansion type
  typedef std::vector<real> multipole_type;
  //! Local expansion type
  typedef std::vector<real> local_type;

  /** Initialize a multipole expansion with the size of a box and level number
   * (Optional: If not implemented, default constructor for multipole_type used)
   *
   * @param[in] M The multipole to be initialized
   * @param[in] extents The dimensions of the box containing the multipole
   * @param[in] level The level number of the box. 0: Root box
   */
  void init_multipole(multipole_type& M,
                      const point_type& extents, unsigned level) const {
    (void) M;
    (void) extents;
    (void) level;
  }
  /** Initialize a local expansion with the size of a box at this level
   * (Optional: If not implemented, default constructor for local_type used)
   *
   * @param[in] L The local expansion to be initialized
   * @param[in] extents The dimensions of the box containing the expansion
   * @param[in] level The level number of the box. 0: Root box
   */
  void init_local(local_type& L,
                  const point_type& extents, unsigned level) const {
    (void) L;
    (void) extents;
    (void) level;
  }

  /** Kernel evaluation
   * K(t,s) where s is the source
   *          and t is the target.
   *
   * @param[in] t,s The target and source to evaluate the kernel
   * @return The Kernel evaluation, K(t,s)
   */
  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    (void) t;
    (void) s;
    return kernel_value_type(0);
  }

  /** Kernel P2M operation
   * M += c * Op(s) where M is the multipole,
   *                      s is the source,
   *                  and c is the charge.
   *
   * @param[in] source The source to accumulate into the multipole
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    (void) source;
    (void) charge;
    (void) center;
    (void) M;
  }

  /** Kernel M2M operation
   * M_t += Op(M_s) where M_t is the target
   *                  and M_s is the source.
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Msource includes the influence of all sources within its box
   */
  void M2M(const multipole_type& source,
           multipole_type& target,
           const point_type& translation) const {
    (void) source;
    (void) target;
    (void) translation;
  }

  /** Kernel M2P operation
   * r += Op(M, t) where M is the multipole,
   *                     t is the target,
   *                 and r is the result.
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] target The target at which to evaluate the multipole
   * @param[in,out] result The corresponding result value to accumulate into
   * @pre M includes the influence of all sources within its box
   */
  void M2P(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result) const {
    (void) M;
    (void) center;
    (void) target;
    (void) result;
  }

  /** Kernel M2L operation
   * L += Op(M) where L is the local expansion
   *              and M is the multipole expansion
   *
   * @param[in] source The multpole expansion source
   * @param[in,out] target The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre source includes the influence of all sources within its box
   */
  void M2L(const multipole_type& source,
           local_type& target,
           const point_type& translation) const {
    (void) source;
    (void) target;
    (void) translation;
  }

  /** Kernel L2L operation
   * L_t += Op(L_s) where L_t is the target
   *                  and L_s is the source.
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre source includes the influence of all sources outside its box
   */
  void L2L(const local_type& source,
           local_type& target,
           const point_type& translation) const {
    (void) source;
    (void) target;
    (void) translation;
  }

  /** Kernel L2P operation
   * r += Op(L, t) where L is the local expansion,
   *                     t is the target,
   *                 and r is the result.
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] target The target of this L2P operation
   * @param[in] result The corresponding result value to accumulate into
   * @pre L includes the influence of all sources outside its box
   */
  void L2P(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
    (void) L;
    (void) center;
    (void) target;
    (void) result;
  }


#if 0
  /*******************************************************************/
  /************************* Optional ********************************/
  /*******************************************************************/

  /* The methods below may be implemented to potentially optimize the P2P
   * operations. If these methods are not implemented, the P2P will be delegated
   * to the Direct.hpp methods which use K.operator()(target_type,source_type)
   * for Kernel evaluations.
   */

  /** Optional Kernel value source and target transposition
   * K(t,s) -> K(s,t)
   * Often, a kernel has a symmetry in s and t that can be computed faster than
   * by calling the evaluation operator. If this function is implemented, the
   * computation may use it to prevent uneccessary calls to the evaluation
   * operator and accelerate the P2P interactions.
   *
   * @param[in] kst A kernel value that was returned from operator()(t,s)
   * @returns The value of operator()(s,t)
   */
  kernel_value_type transpose(const kernel_value_type& kst) const {
    return kst;
  }

  /** Optional Kernel vectorized non-symmetric P2P operation
   * r_i += sum_j K(t_i,s_j) * c_j
   *
   * @param[in] s_first,s_last Iterator pair to the sources
   * @param[in] c_first Iterator to charges corresponding to sources
   * @param[in] t_first,t_last Iterator pair to the targets
   * @param[in] r_first Iterator to result accumulators corresponding to targets
   */
  template <typename SourceIter, typename ChargeIter,
            typename TargetIter, typename ResultIter>
  void P2P(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first) const {
    (void) s_first;
    (void) s_last;
    (void) c_first;
    (void) t_first;
    (void) t_last;
    (void) r_first;
  }

  /** Optional Kernel vectorized symmetric P2P operation
   * This can occur when s_i == t_i for all i.
   * rt_i += sum_j K(t_i, s_j) * cs_j
   * rs_j += sum_i K(s_j, t_i) * ct_i
   *
   * @param[in] s_first,s_last Iterator pair to the sources
   * @param[in] cs_first Iterator to charges corresponding to sources
   * @param[in] t_first,t_last Iterator pair to the targets
   * @param[in] ct_first Iterator to charges corresponding to targets
   * @param[in] rt_first Iterator to result accumulators corresponding to targets
   * @param[in] rs_first Iterator to result accumulators corresponding to sources
   */
  template <typename SourceIter, typename ChargeIter,
            typename TargetIter, typename ResultIter>
  void P2P(SourceIter s_first, SourceIter s_last, ChargeIter cs_first,
           TargetIter t_first, TargetIter t_last, ChargeIter ct_first,
           ResultIter rt_first, ResultIter rs_first) const {
    (void) s_first;
    (void) s_last;
    (void) cs_first;
    (void) t_first;
    (void) t_last;
    (void) ct_first;
    (void) rt_first;
    (void) rs_first;
  }


  /** Optional Kernel vectorized P2M operation
   * M = sum_j Op(s_j) * c_j where M is the multipole and s_j are the sources
   *
   * @param[in] s_first,s_last Iterator pair to the sources
   * @param[in] c_first Iterator to charges corresponding to sources
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   * @pre M is the result of init_multipole
   */
  template <typename SourceIter, typename ChargeIter>
  void P2M(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           const point_type& center, multipole_type& M) const {
    (void) s_first;
    (void) s_last;
    (void) c_first;
    (void) center;
    (void) M;
  }

  /** Optional Kernel vectorized M2P operation
   * r_i += Op(t_i, M) where M is the multipole and r_i are the results
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] t_first,t_last Iterator pair to the targets
   * @param[in] r_first Iterator to the result accumulator
   * @pre M includes the influence of all sources within its box
   */
  template <typename TargetIter, typename ResultIter>
  void M2P(const multipole_type& M, const point_type& center,
           TargetIter t_first, TargetIter t_last,
           ResultIter r_first) const {
    (void) M;
    (void) center;
    (void) t_first;
    (void) t_last;
    (void) r_first;
  }

  /** Optional Kernel vectorized L2P operation
   * r_i += Op(t_i, L) where L is the local expansion and r_i are the results
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] t_first,t_last Iterator pair to the targets
   * @param[in] r_first Iterator to the result accumulator
   * @pre L includes the influence of all sources outside its box
   */
  template <typename TargetIter, typename ResultIter>
  void L2P(const local_type& L, const point_type& center,
           TargetIter t_first, TargetIter t_last,
           ResultIter r_first) const {
    (void) L;
    (void) center;
    (void) t_first;
    (void) t_last;
    (void) r_first;
  }
#endif  // By default, exclude optional methods
};
