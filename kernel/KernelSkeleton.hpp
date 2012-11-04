#pragma once
/** @file KernelSkeleton
 * @brief An example Kernel implementation that explains the required and
 * optionals methods and types to be included in a Kernel class.
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
  void init_multipole(multipole_type& M, double box_size) {
    (void) M;
    (void) box_size;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, double box_size) {
    (void) L;
    (void) box_size;
  }

  /** Kernel evaluation
   * K(t,s) where s is the source point and t is the target point
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
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The point source
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void P2M(const point_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    (void) source;
    (void) charge;
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

  /** Kernel M2P operation
   * r += Op(M) where M is the multipole and r is the result
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] target The target point position
   * @param[in,out] result The target's corresponding result to accumulate into
   * @pre M includes the influence of all points within its box
   */
  void M2P(const multipole_type& M, const point_type& center,
           const point_type& target, result_type& result) const {
    (void) M;
    (void) center;
    (void) target;
    (void) result;
  }

  /** Kernel M2L operation
   * L += Op(M) where L is the local expansion and M is the multipole
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
   * r += Op(L) where L is the local expansion and r is the result
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre L includes the influence of all points outside its box
   */
  void L2P(const local_type& L, const point_type& center,
           const point_type& target, result_type& result) const {
    (void) L;
    (void) center;
    (void) target;
    (void) result;
  }


  /*******************************************************************/
  /************************* Optional ********************************/
  /*******************************************************************/
  /* The methods below may be implemented to potentially optimize the P2P
   * operations. If these methods are not implemented, the P2P will be delegated
   * to the Direct.hpp methods which use K.operator()(point_type,point_type) for
   * Kernel evaluations.
   */


  /** Optional Kernel value source and target transposition
   * K(t,s) -> K(s,t)
   * Often, a kernel has a symmetry in s and t that can be computed faster than
   * by calling the evaluation operator. If this function is implemented, the
   * computation may use it to prevent uneccessary calls to the evaluation
   * operator and accelerate the P2P evaluations and
   *
   * @param[in] kst A kernel value that was returned from operator()(s,t)
   * @returns The value of K(t,s)
   */
  kernel_value_type transpose(const kernel_value_type& kst) const {
    return kst;
  }

  /** Optional Kernel vectorized non-symmetric P2P operation
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

  /** Optional Kernel vectorized symmetric P2P operation
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


  /** Optional Kernel vectorized P2M operation
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

  /** Optional Kernel vectorized M2P operation
   * r_i += Op(M) where M is the multipole and r_i are the results
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

  /** Optional Kernel vectorized L2P operation
   * r_i += Op(L) where L is the local expansion and r_i are the results
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
};
