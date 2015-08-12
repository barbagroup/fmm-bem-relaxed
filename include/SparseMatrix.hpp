#pragma once

/**
 * Simple sparse matrix class (CSR)
 */

#include <vector>

template <typename I, typename T>
struct SparseMatrix
{
  typedef T value_type;
  typedef I index_type;
  // data
  I rows, cols, nnz;
  std::vector<I> offsets, indices;
  std::vector<T> vals;

  // default constructor
  SparseMatrix() : rows(0), cols(0), nnz(0), offsets(0), indices(0), vals(0) {};
  // empty matrix
  SparseMatrix(int r, int c, int nz)
    : rows(r), cols(c), nnz(nz), offsets(r+1), indices(nz), vals(nz) {};

  // matvec
  template <typename VecType>
  std::vector<T> dot(VecType& x) const
  {
    // init return vector
    std::vector<T> y(x.size(),T(0));

    // matvec
    I jj, j;
    T yy;
    #pragma omp parallel for private(j,jj,yy)
    for (I i=0; i<rows; i++) {
      yy = T(0);
      for (jj=offsets[i]; jj<offsets[i+1]; jj++) {
        j = indices[jj];

        // accumulate into register
        yy += x[j]*vals[jj];
      }
      // write out
      y[i] = yy;
    }
    return y;
  }

  // matvec with droptol
  template <typename VecType>
  std::vector<T> dot(VecType& x, double droptol) const
  {
    // init return vector
    std::vector<T> y(x.size(),T(0));

    // matvec
    I j, jj;
    T yy, v;
    #pragma omp parallel for private(j,jj,yy,v)
    for (I i=0; i<rows; i++) {
      yy = T(0);
      for (jj=offsets[i]; jj<offsets[i+1]; jj++) {
        j = indices[jj];

        v = vals[jj];
        // accumulate large enough values into register
        yy += x[j]*(fabs(v) >= droptol)*v;
      }
      // write out
      y[i] = yy;
    }
    return y;
  }
  // assignment operator
  SparseMatrix& operator=(const SparseMatrix& m) {
    this->rows = m.rows;
    this->cols = m.cols;
    this->nnz = m.nnz;

    this->offsets.resize(m.rows+1);
    this->indices.resize(m.nnz);
    this->vals.resize(m.nnz);
    this->offsets = m.offsets;
    this->indices = m.indices;
    this->vals = m.vals;

    return *this;
  }

  void resize(int r, int c, int nz) {
    rows = r;
    cols = c;
    nnz = nz;
    offsets.resize(r+1);
    indices.resize(nnz);
    vals.resize(nnz);
  }
  // copy constructor
  SparseMatrix(const SparseMatrix& m) : rows(m.rows), cols(m.cols), nnz(m.nnz), offsets(m.offsets), indices(m.nnz), vals(m.vals) {};

  // destructor
  ~SparseMatrix() {
    offsets.resize(0);
    indices.resize(0);
    vals.resize(0);
  }

  // return storage size in bytes
  auto storage_size() -> decltype(sizeof(T) + sizeof(I))
  {
    auto index_storage = (rows+1 + nnz + 3)*sizeof(I);
    auto value_storage = nnz * sizeof(T);
    return index_storage + value_storage;
  }
};

template <typename Matrix, typename Vector>
std::vector<typename Matrix::value_type> matvec(const Matrix& A, const Vector& x)
{
  return A.dot(x);
}

template <typename Matrix, typename Vector>
std::vector<typename Matrix::value_type> matvec(const Matrix& A, const Vector& x, const double droptol)
{
  return A.dot(x,droptol);
}

