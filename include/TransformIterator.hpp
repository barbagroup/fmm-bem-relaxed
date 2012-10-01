#pragma once

#include <type_traits>

template <class UnaryFunction,
          class Iterator>
class transform_iterator
{
public:
  typedef std::result_of<const UnaryFunction(iterator_traits<Iterator>::reference)>::type reference;
  typedef std::remove_cv<std::remove_reference<reference>>::type value_type;
  typedef std::add_pointer<value_type>::type pointer;

  typedef iterator_traits<Iterator>::difference_type difference_type;
  typedef iterator_traits<Iterator>::iterator_category iterator_category;

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
    return ++it_;
  }
  transform_iterator& operator--() {
    return --it_;
  }

private:
  Iterator it_;
  UnaryFunction f_;
};
