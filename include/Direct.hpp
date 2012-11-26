#pragma once
/** @file Direct.hpp
 * @brief Dispatch methods for P2P stage
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>
#include <iterator>
#include <vector>

class Direct
{
  /** If no other Direct dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    std::cerr << "Kernel does not have a correct op() or P2P!\n";
    exit(1);
  }

  /** Non-symmetric P2P dispatch
   */
  template <typename Kernel,
	    typename SourceIter, typename ChargeIter,
	    typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_asymm>::type
  eval(const Kernel& K,
	 SourceIter s_begin, SourceIter s_end,
	 ChargeIter c_begin,
	 TargetIter t_begin, TargetIter t_end,
	 ResultIter r_begin)
  {
    K.P2P(s_begin, s_end, c_begin,
	  t_begin, t_end, r_begin);
  }

  /** Non-symmetric P2P using the evaluation operator
   * r_i += sum_j K(t_i, s_j) * c_j
   *
   * @param[in] ...
   */
  template <typename Kernel,
	    typename SourceIter, typename ChargeIter,
	    typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &&
			  !KernelTraits<Kernel>::has_vector_P2P_asymm>::type
  eval(const Kernel& K,
       SourceIter s_begin, SourceIter s_end,
       ChargeIter c_begin,
       TargetIter t_begin, TargetIter t_end,
       ResultIter r_begin)
  {
    // TODO?
    // Optimize on if(std::iterator_traits<All Iters>::iterator_category == random_access_iterator)
    // to eliminate multiple increments

    typedef typename std::iterator_traits<ResultIter>::reference result_reference;

    for ( ; t_begin!=t_end; ++t_begin, ++r_begin) {
      result_reference r = *r_begin;

      SourceIter s = s_begin;
      ChargeIter c = c_begin;
      for ( ; s != s_end; ++s, ++c)
        r += K(*t_begin, *s) * (*c);
    }
  }

  /** Symmetric P2P dispatch
   */
  template <typename Kernel,
	    typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
	 SourceIter p1_begin, SourceIter p1_end,
	 ChargeIter c1_begin,
	 SourceIter p2_begin, SourceIter p2_end,
	 ChargeIter c2_begin,
	 ResultIter r1_begin, ResultIter r2_begin)
  {
    K.P2P(p1_begin, p1_end, c1_begin,
	  p2_begin, p2_end, c2_begin,
	  r1_begin, r2_begin);
  }


  /** Symmetric P2P with only kernel evaluations
   * r2_i += sum_j K(p2_i, p1_j) * c1_j
   * r1_j += sum_i K(p1_j, p2_i) * c2_i
   *
   * @param[in] ...
   */
  template <typename Kernel,
	    typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &&
			  !KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
       SourceIter p1_begin, SourceIter p1_end,
       ChargeIter c1_begin,
       SourceIter p2_begin, SourceIter p2_end,
       ChargeIter c2_begin,
       ResultIter r1_begin, ResultIter r2_begin)
  {
    typedef typename std::iterator_traits<ResultIter>::reference result_reference;

    // TODO
    // Optimize on p1 == p2 (i.e. self-box interaction)
    // Optimize on random_access_iterator?

    for ( ; p1_begin != p1_end; ++p1_begin, ++c1_begin, ++r1_begin) {
      result_reference r1 = *r1_begin;

      SourceIter p2i = p2_begin;
      ChargeIter c2i = c2_begin;
      ResultIter r2i = r2_begin;
      for ( ; p2i != p2_end; ++p2i, ++c2i, ++r2i) {
        r1   += K(*p1_begin, *p2i) * (*c2i);
        *r2i += K(*p2i, *p1_begin) * (*c1_begin);
      }
    }
  }

  /** Symmetric P2P with transpose
   * r2_i += sum_j K(p2_i, p1_j) * c1_j
   * r1_j += sum_i K(p1_j, p2_i) * c2_i
   *
   * @param[in] ...
   */
  template <typename Kernel,
	    typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &&
			  KernelTraits<Kernel>::has_transpose &&
			  !KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
	 SourceIter p1_begin, SourceIter p1_end,
	 ChargeIter c1_begin,
	 SourceIter p2_begin, SourceIter p2_end,
	 ChargeIter c2_begin,
	 ResultIter r1_begin, ResultIter r2_begin)
  {
    typedef typename std::iterator_traits<ResultIter>::reference result_reference;
    typedef typename Kernel::kernel_value_type kernel_value_type;

    // TODO
    // Optimize on p1 == p2 (i.e. self-box interaction)
    // Optimize on random_access_iterator?

    for ( ; p1_begin != p1_end; ++p1_begin, ++c1_begin, ++r1_begin) {
      result_reference r1 = *r1_begin;

      SourceIter p2i = p2_begin;
      ChargeIter c2i = c2_begin;
      ResultIter r2i = r2_begin;
      for ( ; p2i != p2_end; ++p2i, ++c2i, ++r2i) {
        kernel_value_type k12 = K(*p1_begin, *p2i);
        r1   += k12 * (*c2i);
        *r2i += K.transpose(k12) * (*c1_begin);
      }
    }
  }

 public:

  /** Non-symmetric P2P dispatch
   */
  template <typename Kernel,
	    typename SourceIter, typename ChargeIter,
	    typename TargetIter, typename ResultIter>
  inline static void matvec(const Kernel& K,
                            SourceIter s_begin, SourceIter s_end,
			    ChargeIter c_begin,
                            TargetIter t_begin, TargetIter t_end,
			    ResultIter r_begin)
  {
    Direct::eval(K,
		 s_begin, s_end, c_begin,
		 t_begin, t_end, r_begin);
  }

  /** Symmetric P2P dispatch
   */
  template <typename Kernel,
	    typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(const Kernel& K,
                            SourceIter p1_begin, SourceIter p1_end,
			    ChargeIter c1_begin,
                            SourceIter p2_begin, SourceIter p2_end,
			    ChargeIter c2_begin,
                            ResultIter r1_begin, ResultIter r2_begin)
  {
    Direct::eval(K,
		 p1_begin, p1_end, c1_begin,
		 p2_begin, p2_end, c2_begin,
		 r1_begin, r2_begin);
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
    Direct::matvec(K,
                   s.begin(), s.end(), c.begin(),
                   t.begin(), t.end(), r.begin());
  }

  /** Convenience function for std::vector when s == t == p
   */
  template <typename Kernel>
  inline static void matvec(const Kernel& K,
                            const std::vector<typename Kernel::source_type>& p,
                            const std::vector<typename Kernel::charge_type>& c,
                            std::vector<typename Kernel::result_type>& r)
  {
    // TODO
    Direct::matvec(K, p, c, p, r);
  }
};
