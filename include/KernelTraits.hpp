#pragma once

#include <iterator>
#include <iostream>

// Automatically derive !=, <=, >, and >= from a class's == and <
#include <utility>
using namespace std::rel_ops;

#define SFINAE_TEMPLATE(NAME, OP)                                       \
  template <typename ReturnType, typename... Args>                      \
  struct NAME {                                                         \
    template <class A, ReturnType (A::*)(Args...) const> struct SFINAE {}; \
    template <class A> static std::true_type  sfinae(SFINAE<A,&A::OP>*); \
    template <class A> static std::false_type sfinae(...);              \
    static constexpr bool value = decltype(sfinae<Kernel>(0))::value;   \
  }


template <typename Kernel>
struct KernelTraits {
  typedef KernelTraits<Kernel> self_type;

  typedef Kernel kernel_type;

  typedef typename kernel_type::source_type       source_type;
  typedef typename kernel_type::target_type       target_type;
  typedef typename kernel_type::charge_type       charge_type;
  typedef typename kernel_type::kernel_value_type kernel_value_type;
  typedef typename kernel_type::result_type       result_type;

  // A dummy iterator adaptor to check for vectorized methods
  template <typename T>
  struct dummy_iterator {
    typedef T                          value_type;
    typedef T&                         reference;
    typedef T*                         pointer;
    typedef std::ptrdiff_t             difference_type;
    typedef std::forward_iterator_tag  iterator_category;

    dummy_iterator() {}
    dummy_iterator(const dummy_iterator&) {}
    dummy_iterator& operator=(const dummy_iterator&) { return *this; }

    bool operator==(const dummy_iterator&) const { return true; }
    dummy_iterator& operator++() { return *this; }
    value_type& operator*() { return value; }
    const value_type& operator*() const { return value; }
    pointer operator->() { return nullptr; }
    value_type value;
  };

  typedef dummy_iterator<source_type> source_iterator;
  typedef dummy_iterator<charge_type> charge_iterator;
  typedef dummy_iterator<target_type> target_iterator;
  typedef dummy_iterator<result_type> result_iterator;

  // Kernel Evaluations and P2P
  SFINAE_TEMPLATE(HasEvalOp,operator());
  static constexpr bool has_eval_op =
      HasEvalOp<kernel_value_type,
                const target_type&, const source_type&>::value;
  SFINAE_TEMPLATE(HasTranspose,transpose);
  static constexpr bool has_transpose =
      HasTranspose<kernel_value_type,
                   const kernel_value_type&>::value;
  SFINAE_TEMPLATE(HasP2P,P2P);
  static constexpr bool has_vector_P2P_symm =
      HasP2P<void,
             source_iterator, source_iterator, charge_iterator,
             target_iterator, target_iterator, charge_iterator,
             result_iterator, result_iterator>::value;
  static constexpr bool has_vector_P2P_asymm =
      HasP2P<void,
             source_iterator, source_iterator, charge_iterator,
             target_iterator, target_iterator, result_iterator>::value;

  static constexpr bool is_valid_kernel = has_eval_op || has_vector_P2P_asymm;

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << "has_eval_op: " << traits.has_eval_op << std::endl;
    s << "has_transpose: " << traits.has_transpose << std::endl;
    s << "has_vector_P2P_symm: " << traits.has_vector_P2P_symm << std::endl;
    s << "has_vector_P2P_asymm: " << traits.has_vector_P2P_asymm << std::endl;
    return s;
  }
};




template <typename Kernel>
struct ExpansionTraits : public KernelTraits<Kernel>
{
  typedef ExpansionTraits<Kernel> self_type;
  typedef KernelTraits<Kernel> super_type;

  typedef typename super_type::source_type       source_type;
  typedef typename super_type::target_type       target_type;
  typedef typename super_type::charge_type       charge_type;
  typedef typename super_type::kernel_value_type kernel_value_type;
  typedef typename super_type::result_type       result_type;

  typedef typename super_type::source_iterator source_iterator;
  typedef typename super_type::charge_iterator charge_iterator;
  typedef typename super_type::target_iterator target_iterator;
  typedef typename super_type::result_iterator result_iterator;


  static constexpr unsigned dimension = Kernel::dimension;
  typedef typename Kernel::point_type        point_type;

  // TODO: Check that dimension and point_type make sense

  typedef typename Kernel::multipole_type    multipole_type;
  typedef typename Kernel::local_type        local_type;

  // Initializers
  SFINAE_TEMPLATE(HasInitMultipole,init_multipole);
  static constexpr bool has_init_multipole =
      HasInitMultipole<void,
                       multipole_type&, const point_type&, unsigned>::value;
  SFINAE_TEMPLATE(HasInitLocal,init_local);
  static constexpr bool has_init_local =
      HasInitLocal<void,
                   local_type&, const point_type&, unsigned>::value;

  // Kernel Evaluations and P2P
  static constexpr bool has_eval_op          = super_type::has_eval_op;
  static constexpr bool has_transpose        = super_type::has_transpose;
  static constexpr bool has_vector_P2P_symm  = super_type::has_vector_P2P_symm;
  static constexpr bool has_vector_P2P_asymm = super_type::has_vector_P2P_asymm;

  // P2M
  SFINAE_TEMPLATE(HasP2M,P2M);
  static constexpr bool has_P2M =
      HasP2M<void,
             const source_type&, const charge_type&,
             const point_type&, multipole_type&>::value;
  static constexpr bool has_vector_P2M =
      HasP2M<void,
             source_iterator, source_iterator, charge_iterator,
             const point_type&, multipole_type&>::value;

  // M2M
  SFINAE_TEMPLATE(HasM2M,M2M);
  static constexpr bool has_M2M =
      HasM2M<void,
             const multipole_type&, multipole_type&, const point_type&>::value;

  // M2P
  SFINAE_TEMPLATE(HasM2P,M2P);
  static constexpr bool has_M2P =
      HasM2P<void,
             const multipole_type&, const point_type&,
             const target_type&, result_type&>::value;
  static constexpr bool has_vector_M2P =
      HasM2P<void,
             const multipole_type&, const point_type&,
             target_iterator, target_iterator, result_iterator>::value;

  // M2L
  SFINAE_TEMPLATE(HasM2L,M2L);
  static constexpr bool has_M2L =
      HasM2L<void,
             const multipole_type&, local_type&, const point_type&>::value;

  // L2L
  SFINAE_TEMPLATE(HasL2L,L2L);
  static constexpr bool has_L2L =
      HasL2L<void,
             const local_type&, local_type&, const point_type&>::value;

  // L2P
  SFINAE_TEMPLATE(HasL2P,L2P);
  static constexpr bool has_L2P =
      HasL2P<void,
             const local_type&, const point_type&,
             const target_type&, result_type&>::value;
  static constexpr bool has_vector_L2P =
      HasL2P<void,
             const local_type&, const point_type&,
             target_iterator, target_iterator, result_iterator>::value;

  static constexpr bool is_valid_treecode =
      (has_eval_op || has_vector_P2P_asymm) &&
      (has_P2M || has_vector_P2M) &&
      (has_M2M) &&
      (has_M2P || has_vector_M2P);

  static constexpr bool is_valid_fmm = (has_eval_op || has_vector_P2P_asymm) &&
      (has_P2M || has_vector_P2M) &&
      (has_M2M) &&
      (has_M2L) &&
      (has_L2L) &&
      (has_L2P || has_vector_L2P);

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << static_cast<super_type>(traits);
    s << "has_init_multipole: " << traits.has_init_multipole << std::endl;
    s << "has_init_local: " << traits.has_init_local << std::endl;
    s << "has_vector_P2M: " << traits.has_vector_P2M << std::endl;
    s << "has_P2M: " << traits.has_P2M << std::endl;
    s << "has_vector_P2M: " << traits.has_vector_P2M << std::endl;
    s << "has_M2M: " << traits.has_M2M << std::endl;
    s << "has_M2P: " << traits.has_M2P << std::endl;
    s << "has_vector_M2P: " << traits.has_vector_M2P << std::endl;
    s << "has_M2L: " << traits.has_M2L << std::endl;
    s << "has_L2L: " << traits.has_L2L << std::endl;
    s << "has_L2P: " << traits.has_L2P << std::endl;
    s << "has_vector_L2P: " << traits.has_vector_L2P << std::endl;
    return s;
  }
};

#undef SFINAE_TEMPLATE
