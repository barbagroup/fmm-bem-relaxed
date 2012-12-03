#pragma once

template <typename Tree, typename Indexable>
struct body_map {
  typedef body_map<Indexable> self_type;
  typedef typename Indexable::value_type value_type;
  typedef typename Tree::box_type box_type;
  typedef typename Tree::body_type body_type;
  typedef typename Tree::body_iterator body_iterator;

  // Constructor
  body_map(const Indexable& body_value)
    : body_value_(body_value) {
  }
  // RValue constructor for efficiency
  body_map(Indexable&& body_value)
    : body_value_(std::move(body_value)) {
  }

  const value_type& operator()(const body_type& b) const {
    return body_value_[b.number()];  // TODO TEMP: number to workaround permute
  }
  value_type& operator()(const body_type& b) {
    return box_value_[b.number()];  // TODO TEMP: number to workaround permute
  }

  // self_type& to prevent copying, Indexable may be a heavy container
  typedef transform_iterator<body_iterator, self_type&> body_value_iterator;

  body_value_iterator value_begin(const box_type& box) const {
    return make_transform_iterator<self_type&>(b.body_begin(), *this);
  }
  body_value_iterator value_end(const box_type& box) const {
    return make_transform_iterator<self_type&>(b.body_end(), *this);
  }

private:
  Indexable body_value_;
};


template <typename Tree, typename Indexable>
struct box_map {
  typedef box_map<Indexable> self_type;
  typedef typename Indexable::value_type value_type;
  typedef typename Tree::box_type box_type;

  // Constructor
  box_map(const Indexable& body_value)
    : body_value_(body_value) {
  }
  // RValue constructor for efficiency
  box_map(Indexable&& body_value)
    : body_value_(std::move(body_value)) {
  }

  const value_type& operator()(const box_type& b) const {
    return box_value_[b.index()];
  }
  value_type& operator()(const box_type& b) {
    return box_value_[b.index()];
  }

private:
  Indexable box_value_;
};
