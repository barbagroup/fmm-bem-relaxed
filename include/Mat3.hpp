#pragma once

/**
 * Simple 3x3 matrix class for analytical / semi-analytical integrals
 */

template <typename T>
struct Mat3 {
  // 3 rows, 3 cols, 9 values
  T vals_[9];

  Mat3() { for (unsigned i=0; i<9; i++) vals_[i] = 0.; };
  template <typename IterType>
  Mat3(IterType start, IterType end) {
    assert(end-start == 9); // ensure correct # values
    unsigned i=0;
    for ( ; start!=end; ++start, ++i) vals_[i] = *start;
  };
  // copy constructor
  Mat3(const Mat3<T>& M) {
    for (unsigned i=0; i<9; i++) vals_[i] = M.vals_[i];
  }
  // initialise with a uniform value
  Mat3(double v) {
    for (unsigned i=0; i<9; i++) vals_[i] = v;
  }

  // return negated matrix
  Mat3<T> operator-() const {
    Mat3<T> temp;
    for (unsigned i=0; i<9; i++) temp.vals_[i] = -vals_[i];
    return temp;
  };
  // accessors
  const T& operator()(unsigned i, unsigned j) const {
    return vals_[i*3+j];
  };
  T& operator()(unsigned i, unsigned j) {
    return vals_[i*3+j];
  };
  // M = M2
  Mat3<T>& operator=(const Mat3<T>& M) {
    if (this == &M) return *this;
    for (unsigned i=0; i<9; i++) this->vals_[i] = M.vals_[i];
    return *this;
  }
  // M += M2
  Mat3<T>& operator+=(const Mat3<T>& M) {
    for (unsigned i=0; i<9; i++) this->vals_[i] += M.vals_[i];

    return *this;
  }
  // M = M1 + M2
  Mat3<T> operator+(const Mat3<T>& M) {
    Mat3<T> temp(*this);

    return (temp += M);
  }
  Vec<3,T> operator*(const Vec<3,T>& x) const {
    Vec<3,T> result;
    result[0] = vals_[0]*x[0]+vals_[1]*x[1]+vals_[2]*x[2];
    result[1] = vals_[3]*x[0]+vals_[4]*x[1]+vals_[5]*x[2];
    result[2] = vals_[6]*x[0]+vals_[7]*x[1]+vals_[8]*x[2];
    return result;
  }
  Mat3<T> operator*(const double x) const {
    Mat3<T> result;

    for (unsigned i=0; i<9; i++) result.vals_[i] = x*vals_[i];

    return result;
  }
  // matvec
  Vec<3,T> multiply(const Vec<3,T>& x) const {
    Vec<3,T> result;
    result[0] = vals_[0]*x[0]+vals_[1]*x[1]+vals_[2]*x[2];
    result[1] = vals_[3]*x[0]+vals_[4]*x[1]+vals_[5]*x[2];
    result[2] = vals_[6]*x[0]+vals_[7]*x[1]+vals_[8]*x[2];
    return result;
  }
  Mat3<T> multiply(const double x) const {
    Mat3<T> M(*this);

    for (unsigned i=0; i<9; i++) M.vals_[i]*=x;
    return M;
  }
};
