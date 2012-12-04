#pragma once

#include <iterator>
#include <type_traits>
#include <utility>

// Automatically derive !=, <=, >, and >= from a class's == and <
using namespace std::rel_ops;


template <class Iterator,
          class UnaryFunction>
class transform_iterator
{
  typedef typename std::iterator_traits<Iterator>::reference it_reference;
public:
  typedef typename std::result_of<const UnaryFunction(it_reference)>::type reference;
  typedef typename std::remove_cv<std::remove_reference<reference>>::type value_type;
  typedef typename std::add_pointer<value_type>::type pointer;

  typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;

  // CONSTRUCTORS

  // Construct an invalid iterator
  transform_iterator() {}
  // Construct a valid iterator
  transform_iterator(const Iterator& x, const UnaryFunction& f)
      : it_(x), f_(f) {
  }
  // Copy constructor
  transform_iterator(const transform_iterator& t)
      : it_(t.it_), f_(t.f_) {
  }

  // ACCESSORS

  const UnaryFunction& functor() const {
    return f_;
  }
  const Iterator& base() const {
    return it_;
  }

  // OPERATORS

  const reference operator*() const {
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


template <class F, typename IT>
transform_iterator<IT,F> make_transform_iterator(IT it, F f) {
  return transform_iterator<IT,F>(it, f);
}
