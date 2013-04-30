#pragma once
/** @file Direct.hpp
 * @brief Dispatch methods for P2P stage
 *
 */

#include "KernelTraits.hpp"

#include <type_traits>
#include <iterator>
#include <vector>

class Direct {
	/** If no other Direct dispatcher matches */
	template <typename... Args>
	inline static void eval(Args...) {
		std::cerr << "Kernel does not have a correct op() or P2P!\n";
		exit(1);
	}

	/** Asymmetric P2P dispatch
	 */
	template <typename Kernel,
	          typename SourceIter, typename ChargeIter,
	          typename TargetIter, typename ResultIter>
	inline static
	typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_asymm>::type
	eval(const Kernel& K,
	     SourceIter s_first, SourceIter s_last, ChargeIter c_first,
	     TargetIter t_first, TargetIter t_last, ResultIter r_first)
	{
		K.P2P(s_first, s_last, c_first,
		      t_first, t_last, r_first);
	}

  /** Symmetric P2P, off-diagonal block dispatch
   * @pre source_type == target_type
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
       SourceIter p1_first, SourceIter p1_last, ChargeIter c1_first,
       ResultIter r1_first,
       SourceIter p2_first, SourceIter p2_last, ChargeIter c2_first,
       ResultIter r2_first)
  {
    K.P2P(p1_first, p1_last, c1_first,
          p2_first, p2_last, c2_first,
          r1_first, r2_first);
  }

  /** Dual-Evaluation dispatch
   */
  template <typename Kernel>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &
                          !KernelTraits<Kernel>::has_transpose>::type
  eval(const Kernel& K,
       const typename Kernel::source_type& p1,
       const typename Kernel::charge_type& c1,
       typename Kernel::result_type& r1,
       const typename Kernel::source_type& p2,
       const typename Kernel::charge_type& c2,
       typename Kernel::result_type& r2)
  {
    r1 += K(p1,p2) * c2;
    r2 += K(p2,p1) * c1;
  }

  /** Dual-Evaluation dispatch
   */
  template <typename Kernel>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &
                          KernelTraits<Kernel>::has_transpose>::type
  eval(const Kernel& K,
       const typename Kernel::source_type& p1,
       const typename Kernel::charge_type& c1,
       typename Kernel::result_type& r1,
       const typename Kernel::source_type& p2,
       const typename Kernel::charge_type& c2,
       typename Kernel::result_type& r2)
  {
    typedef typename Kernel::kernel_value_type kernel_value_type;

    kernel_value_type k12 = K(p1,p2);
    r1 += k12 * c2;
    kernel_value_type k21 = K.transpose(k12);
    r2 += k21 * c1;
  }

	/** Asymmetric P2P using the evaluation operator
	 * r_i += sum_j K(t_i, s_j) * c_j
	 *
	 * @param[in] ...
	 */
	template <typename Kernel,
	          typename SourceIter, typename ChargeIter,
	          typename TargetIter, typename ResultIter>
	inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &
                          !KernelTraits<Kernel>::has_vector_P2P_asymm>::type
  eval(const Kernel& K,
       SourceIter s_first, SourceIter s_last, ChargeIter c_first,
       TargetIter t_first, TargetIter t_last, ResultIter r_first)
  {
    // TODO?
    // Optimize on if(std::iterator_traits<All Iters>::iterator_category == random_access_iterator)
    // to eliminate multiple increments

    typedef typename Kernel::result_type result_type;
    typedef typename Kernel::target_type target_type;

    for ( ; t_first != t_last; ++t_first, ++r_first) {
      const target_type& t = *t_first;
      result_type& r       = *r_first;

      SourceIter s = s_first;
      ChargeIter c = c_first;
      for ( ; s != s_last; ++s, ++c)
        r += K(t, *s) * (*c);
    }
  }


  /** Symmetric P2P, off-diagonal block
   * r2_i += sum_j K(p2_i, p1_j) * c1_j
   * r1_j += sum_i K(p1_j, p2_i) * c2_i
   *
   * @param[in] ...
   * @pre source_type == target_type
   * @pre For all i,j we have p1_i != p2_j
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<!KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
       SourceIter p1_first, SourceIter p1_last, ChargeIter c1_first,
       ResultIter r1_first,
       SourceIter p2_first, SourceIter p2_last, ChargeIter c2_first,
       ResultIter r2_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    static_assert(std::is_same<source_type,target_type>::value,
                  "source_type != target_type in symmetric P2P");

    // TODO
    // Optimize on random_access_iterator?

    for ( ; p1_first != p1_last; ++p1_first, ++c1_first, ++r1_first) {
      const source_type& p1 = *p1_first;
      const charge_type& c1 = *c1_first;
      result_type& r1       = *r1_first;

      SourceIter p2i = p2_first;
      ChargeIter c2i = c2_first;
      ResultIter r2i = r2_first;
      for ( ; p2i != p2_last; ++p2i, ++c2i, ++r2i)
        Direct::eval(K, p1, c1, r1, *p2i, *c2i, *r2i);
    }
  }

  /** Symmetric P2P, diagonal block
   * r_i += sum_j K(p_i, p_j) * c_j
   *
   * @pre source_type == target_type
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<!KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
       SourceIter p_first, SourceIter p_last,
       ChargeIter c_first, ResultIter r_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    static_assert(std::is_same<source_type, target_type>::value,
                  "source_type != target_type in symmetric P2P");
    static_assert(std::is_same<source_type,
                               typename SourceIter::value_type>::value,
                  "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                               typename ChargeIter::value_type>::value,
                  "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<result_type,
                               typename ResultIter::value_type>::value,
                  "ResultIter::value_type != Kernel::result_type");

    // TODO
    // Optimize on random_access_iterator?

    SourceIter pi = p_first;
    ChargeIter ci = c_first;
    ResultIter ri = r_first;
    for ( ; pi != p_last; ++pi, ++ci, ++ri) {
      const source_type& p = *pi;
      const charge_type& c = *ci;
      result_type& r       = *ri;

      // The off-diagonal elements (use the symmetry)
      SourceIter pj = p_first;
      ChargeIter cj = c_first;
      ResultIter rj = r_first;
      for ( ; pj != pi; ++pj, ++cj, ++rj)
        Direct::eval(K, p, c, r, *pj, *cj, *rj);

      // The diagonal element
      r += K(p, p) * c;
    }
  }



 public:

  /** Asymmetric matvec
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter,
            typename TargetIter, typename ResultIter>
  inline static void matvec(const Kernel& K,
                            SourceIter s_first, SourceIter s_last,
                            ChargeIter c_first,
                            TargetIter t_first, TargetIter t_last,
                            ResultIter r_first)
  {
    Direct::eval(K,
                 s_first, s_last, c_first,
                 t_first, t_last, r_first);
  }

  /** Symmetric matvec, off-diagonal block
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(const Kernel& K,
                            SourceIter p1_first, SourceIter p1_last,
                            ChargeIter c1_first, ResultIter r1_first,
                            SourceIter p2_first, SourceIter p2_last,
                            ChargeIter c2_first, ResultIter r2_first)
  {
    Direct::eval(K,
                 p1_first, p1_last, c1_first, r1_first,
                 p2_first, p2_last, c2_first, r2_first);
  }

  /** Symmetric matvec, diagonal block
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(const Kernel& K,
                            SourceIter p_first, SourceIter p_last,
                            ChargeIter c_first, ResultIter r_first)
  {
    Direct::eval(K,
                 p_first, p_last, c_first, r_first);
  }

  /** Convenience function for std::vector
   */
  template <typename Kernel>
  inline static void matvec(const Kernel& K,
                            const std::vector<typename Kernel::source_type>& s,
                            const std::vector<typename Kernel::charge_type>& c,
                            const std::vector<typename Kernel::target_type>& t,
                            std::vector<typename Kernel::result_type>& r)
  {
    assert(s.size() == c.size());
    assert(t.size() == r.size());

    // Pass to asymmetric matvec
    Direct::matvec(K,
                   s.begin(), s.end(), c.begin(),
                   t.begin(), t.end(), r.begin());
  }

  /** Convenience function for std::vector
   */
  template <typename Kernel>
  inline static void matvec(const Kernel& K,
                            const std::vector<typename Kernel::source_type>& p,
                            const std::vector<typename Kernel::charge_type>& c,
                            std::vector<typename Kernel::result_type>& r)
  {
    assert(p.size() == c.size());
    assert(p.size() == r.size());

    // Pass to symmetric matvec
    Direct::matvec(K,
                   p.begin(), p.end(), c.begin(), r.begin());
  }
};
