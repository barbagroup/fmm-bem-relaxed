#pragma once

/**
 * BLAS and custom matrix operations
 */

/** Very basic matrix type */
template <typename T>
class Matrix
{
 private:
  int rows_, cols_;
  // 1D storage, column-major
  std::vector<T> vals_;
  // T *vals_;

 public:
  Matrix() : rows_(0), cols_(0), vals_(0) {};
  Matrix(int r, int c) : rows_(r), cols_(c) {
    // vals_.resize(rows_*cols_);
    vals_ = std::vector<T>(rows_*cols_);
    // vals_ = new T[rows_*cols_];
  };
  Matrix(const Matrix& M) : rows_(M.rows_), cols_(M.cols_), vals_(M.vals_) {};

  //~Matrix() { delete [] vals_; };

  int rows() { return rows_; };
  int cols() { return cols_; };

  std::vector<T> column(int c) {
    std::vector<T> r(rows_,0);
    for (int i=0; i<rows_; i++) {
      r[i] = vals_[c*rows_+i];
    }
    return r;
    // return std::vector<T>(&vals_[c*rows_],&vals_[c*(rows_+1)]);
  }
  const T& operator()(int i, int j) const {
    return vals_[j*rows_+i];
  }
  T& operator()(int i, int j) {
    return vals_[j*rows_+i];
  }
  // modify parts of the matrix
  void set_column(int col, std::vector<T> v) {
    for (unsigned i=0; i<v.size(); i++) {
      vals_[col*rows_ + i] = v[i];
    }
  }
  // overload for Vec<N,T>
  template <int N, typename T2>
  void set_column(int col, std::vector<Vec<N,T2>> v) {
    for (unsigned i=0; i<v.size(); i++) {
      for (int j=0; j<N; j++) {
        vals_[col*rows_ + N*i + j] = v[i][j];
      }
    }
  }
  void print() {
    for (unsigned i=0; i<rows_; i++) {
      for (unsigned j=0; j<cols_; j++) {
        printf("%.3lg  ",this->operator()(i,j));
      }
      printf("\n");
    }
  }
};

namespace blas {

template <typename T>
T nrm2(std::vector<T>& vec) {
  T res = T(0);
  for (auto it=vec.begin(); it!=vec.end(); ++it) res += (*it) * (*it);
  return std::sqrt(res);
}

template <int N, typename T>
T nrm2(std::vector<Vec<N,T>>& vec) {
  T res = T(0.);

  for (auto it=vec.begin(); it!=vec.end(); ++it) {
    for (int j=0; j<N; j++) {
      res += (*it)[j] * (*it)[j];
    }
  }
  return std::sqrt(res);
}

template <typename T1, typename T2>
void scal(std::vector<T1>& vec, T2 s)
{
  for (auto it=vec.begin(); it!=vec.end(); ++it)
    *it *= s;
}

template <typename T>
T dotc(std::vector<T>& x, std::vector<T>& y) {
  T ret = T(0);
  auto yit = y.begin();
  for (auto it=x.begin(); it!=x.end(); ++it, ++yit)
    ret += (*it) * (*yit);
  return ret;
};

template <typename Matrix, typename T>
void matvec(Matrix& A, std::vector<T>& x, std::vector<T>& r) {
  r = std::vector<T>(x.size(),0);
  unsigned j;
  #pragma omp parallel for private(j)
  for (unsigned i=0; i<A.rows(); i++) {
    for (j=0; j<A.cols(); j++) {
      r[i] += A(i,j)*x[i];
    }
  }
}

/** y <- a*x + y */
template <typename T>
void axpy(std::vector<T> x, std::vector<T>& y, T a) {
  for (unsigned i=0; i<x.size(); i++) {
    y[i] = a*x[i] + y[i];
  }
}

template <int N, typename T>
void axpy(std::vector<Vec<N,T>> x, std::vector<Vec<N,T>>& y, T a) {
  for (unsigned i=0; i<x.size(); i++) {
    y[i] += a*x[i];
  }
}

template <int N, typename T>
void axpy(std::vector<Vec<N,T>> x, std::vector<T>& y, T a) {
  for (unsigned i=0; i<x.size(); i++) {
    for (unsigned j=0; j<N; j++) {
      y[i*N + j] += a*x[i][j];
    }
  }
}

}; // end namespace blas

template <typename Iter>
typename Iter::value_type norm(Iter start, Iter end)
{
  typename Iter::value_type sum = 0;

  for ( ; start != end; ++start) sum += (*start) * (*start);

  return std::sqrt(sum);
}

