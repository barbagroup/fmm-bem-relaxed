#pragma once


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
  void init_multipole(multipole_type& M, double) const {
    M = 0;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, double) const {
    L = 0;
  }

  /** Kernel evaluation
   * K(t,s)
   *
   * @param[in] t,s The target and source points to evaluate the kernel
   */
  kernel_value_type operator()(const point_type& t,
                               const point_type& s) const {
    if (t != s)
      return kernel_value_type(1);
    else
      return kernel_value_type(0);
  }

  /** Kernel P2M operation
   * M = sum_i Op(s_i) * c_i where M is the multipole and s_i are the sources
   *
   * @param[in] p_begin,p_end Iterator pair to the points in this operation
   * @param[in] c_begin Corresponding charge iterator for the points
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   * @pre M is the result of init_multipole
   */
  template <typename PointIter, typename ChargeIter>
  void P2M(PointIter p_begin, PointIter p_end, ChargeIter c_begin,
           const point_type&, multipole_type& M) const {
    for ( ; p_begin != p_end; ++p_begin, ++c_begin)
      M += *c_begin;
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
   * @pre Msource includes the influence of all points within its box
   */
  void M2L(const multipole_type& source,
                 local_type& target,
           const point_type&) const {
    target += source;
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
  void M2P(const multipole_type& M, const point_type&,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      *r_begin += M;
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
           const point_type&) const {
    target += source;
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
  void L2P(const local_type& L, const point_type&,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    for ( ; t_begin!=t_end; ++t_begin, ++r_begin)
      *r_begin += L;
  }
};
