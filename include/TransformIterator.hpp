#pragma once

#include <iterator>
#include <type_traits>
#include <utility>

// Automatically derive !=, <=, >, and >= from a class's == and <
using namespace std::rel_ops;

template <class UnaryFunction,
          class Iterator>
class transform_iterator
{
public:
  typedef typename std::result_of<const UnaryFunction(typename std::iterator_traits<Iterator>::reference)>::type reference;
  typedef typename std::remove_cv<std::remove_reference<reference>>::type value_type;
  typedef typename std::add_pointer<value_type>::type pointer;

  typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;

  // CONSTRUCTORS

  // Construct an invalid iterator
  transform_iterator() {}
  // Construct a valid iterator
  transform_iterator(const Iterator& x, UnaryFunction f)
      : it_(x), f_(f) {
  }

  template<class F2, class I2>
  transform_iterator(const transform_iterator<F2, I2>& t)
      : it_(t.it_), f_(t.f_) {
  }

  // ACCESSORS

  UnaryFunction functor() const {
    return f_;
  }
  const Iterator& base() const {
    return it_;
  }

  // OPERATORS

  reference operator*() const {
    return f_(*it_);
  }
  transform_iterator& operator++() {
    ++it_;
    return *this;
  }
  transform_iterator& operator--() {
    --it_;
    return *this;
  }
  bool operator==(const transform_iterator& it) const {
    return it_ == it.it_;
  }

private:
  Iterator it_;
  UnaryFunction f_;
};


template <typename IT, class F>
transform_iterator<F,IT> make_transform_iter(IT it, F f) {
  return transform_iterator<F,IT>(it, f);
}
