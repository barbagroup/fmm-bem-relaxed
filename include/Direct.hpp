#pragma once
/** @file Direct.hpp
 * @brief Dispatch methods for P2P stage
 *
 */

#include <iterator>
#include <vector>


struct Direct
{
 private:
  /** Type for Kernel::P2P forwarding identification */
  template <bool> struct UseP2P {};
  template <bool> struct UseTranspose {};

  /** Use *Substitution Failure Is Not An Error* and variadic templates
   * to determine if a Kernel has a P2P method with a certain signature
   */
  template <typename K, typename... Args>
  struct HasP2P {
    template <class A, void (A::*)(Args...) const> struct SFINAE {};
    template <class A> static std::true_type  sfinae(SFINAE<A, &A::P2P>*);
    template <class A> static std::false_type sfinae(...);
    static constexpr bool value = decltype(sfinae<K>(0))::value;
  };

  /** Use *Substitution Failure Is Not An Error*
   * to determine if a Kernel has the method
   * K::kernel_value_type K::transpose(const K::kernel_value_type&) const
   */
  template <typename K>
  struct HasTranspose {
    typedef typename K::kernel_value_type kv_type;
    template <class A, kv_type (A::*)(const kv_type&) const> struct SFINAE {};
    template <class A> static std::true_type  sfinae(SFINAE<A, &A::transpose>*);
    template <class A> static std::false_type sfinae(...);
    static constexpr bool value = decltype(sfinae<K>(0))::value;
  };


  /** Forwards any arguments to the Kernel's P2P function
   * @pre HasP2P<Kernel,Args...>::value == true
   */
  template <typename Kernel, typename... Args>
  inline static void matvec(UseP2P<true>, Kernel& K, Args&&... args)
  {
    static_assert(HasP2P<Kernel,Args...>::value,
		  "Attempting to use P2P that doesn't exist!");
    K.P2P(args...);
  }

  /** Non-symmetric P2P
   * r_i += sum_j K(t_i, s_j) * c_j
   *
   * @param[in] ...
   *
   */
  template <typename Kernel,
	    typename PointIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(UseP2P<false>, Kernel& K,
                            PointIter s_begin, PointIter s_end, ChargeIter c_begin,
                            PointIter t_begin, PointIter t_end, ResultIter r_begin)
  {
    // TODO?
    // Optimize on if(std::iterator_traits<All Iters>::iterator_category == random_access_iterator)
    // to eliminate multiple increments

    typedef typename std::iterator_traits<ResultIter>::reference result_reference;

    for ( ; t_begin!=t_end; ++t_begin, ++r_begin) {
      result_reference r = *r_begin;

      PointIter  s = s_begin;
      ChargeIter c = c_begin;
      for ( ; s != s_end; ++s, ++c)
        r += K(*t_begin, *s) * (*c);
    }
  }

  /** Symmetric P2P
   * r2_i += sum_j K(p2_i, p1_j) * c1_j
   * r1_j += sum_i K(p1_j, p2_i) * c2_i
   *
   * @param[in] ...
   */
  template <typename Kernel,
	    typename PointIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(UseP2P<false>, UseTranspose<false>,
                            Kernel& K,
                            PointIter p1_begin, PointIter p1_end, ChargeIter c1_begin,
                            PointIter p2_begin, PointIter p2_end, ChargeIter c2_begin,
                            ResultIter r1_begin, ResultIter r2_begin)
  {
    typedef typename std::iterator_traits<ResultIter>::reference result_reference;

    // TODO
    // Optimize on p1 == p2 (i.e. self-box interaction)
    // Optimize on random_access_iterator?

    for ( ; p1_begin != p1_end; ++p1_begin, ++c1_begin, ++r1_begin) {
      result_reference r1 = *r1_begin;

      auto p2i = p2_begin;
      auto c2i = c2_begin;
      auto r2i = r2_begin;
      for ( ; p2i != p2_end; ++p2i, ++c2i, ++r2i) {
        r1   += K(*p1_begin, *p2i) * (*c2i);
        *r2i += K(*p2i, *p1_begin) * (*c1_begin);
      }
    }
  }

  /** Symmetric P2P
   * r2_i += sum_j K(p2_i, p1_j) * c1_j
   * r1_j += sum_i K(p1_j, p2_i) * c2_i
   *
   * @param[in] ...
   */
  template <typename Kernel,
	    typename PointIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(UseP2P<false>, UseTranspose<true>,
                            Kernel& K,
                            PointIter p1_begin, PointIter p1_end, ChargeIter c1_begin,
                            PointIter p2_begin, PointIter p2_end, ChargeIter c2_begin,
                            ResultIter r1_begin, ResultIter r2_begin)
  {
    typedef typename std::iterator_traits<ResultIter>::reference result_reference;
    typedef typename Kernel::kernel_value_type kernel_value_type;

    // TODO
    // Optimize on p1 == p2 (i.e. self-box interaction)
    // Optimize on random_access_iterator?

    for ( ; p1_begin != p1_end; ++p1_begin, ++c1_begin, ++r1_begin) {
      result_reference r1 = *r1_begin;

      auto p2i = p2_begin;
      auto c2i = c2_begin;
      auto r2i = r2_begin;
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
	    typename PointIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(Kernel& K,
                            PointIter s_begin, PointIter s_end, ChargeIter c_begin,
                            PointIter t_begin, PointIter t_end, ResultIter r_begin)
  {
    typedef HasP2P<Kernel,
                   PointIter, PointIter, ChargeIter,
                   PointIter, PointIter, ResultIter> KernelP2P;

    matvec(UseP2P<KernelP2P::value>(),
           K,
           s_begin, s_end, c_begin,
           t_begin, t_end, r_begin);
  }

  /** Symmetric P2P dispatch
   */
  template <typename Kernel,
	    typename PointIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(Kernel& K,
                            PointIter p1_begin, PointIter p1_end, ChargeIter c1_begin,
                            PointIter p2_begin, PointIter p2_end, ChargeIter c2_begin,
                            ResultIter r1_begin, ResultIter r2_begin)
  {
    typedef HasP2P<Kernel,
                   PointIter, PointIter, ChargeIter,
                   PointIter, PointIter, ChargeIter,
                   ResultIter, ResultIter> KernelP2P;

    typedef HasTranspose<Kernel> KernelTranspose;

    matvec(UseP2P<KernelP2P::value>(), UseTranspose<KernelTranspose::value>(),
           K,
           p1_begin, p1_end, c1_begin,
           p2_begin, p2_end, c2_begin,
           r1_begin, r2_begin);
  }

  /** Convenience function for std::vector
   */
  template <typename Kernel>
  inline static void matvec(Kernel& K,
                            const std::vector<typename Kernel::point_type>& s,
                            const std::vector<typename Kernel::charge_type>& c,
                            const std::vector<typename Kernel::point_type>& t,
                            std::vector<typename Kernel::result_type>& r)
  {
    Direct::matvec(K,
                   s.begin(), s.end(), c.begin(),
                   t.begin(), t.end(), r.begin());
  }

  /** Convenience function for std::vector when s == t == p
   */
  template <typename Kernel>
  inline static void matvec(Kernel& K,
                            const std::vector<typename Kernel::point_type>& p,
                            const std::vector<typename Kernel::charge_type>& c,
                            std::vector<typename Kernel::result_type>& r)
  {
    // TODO
    Direct::matvec(K, p, c, p, r);
  }
};
