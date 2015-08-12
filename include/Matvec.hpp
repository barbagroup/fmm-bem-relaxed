#pragma once

/**
 * matvec for ublas, where:
 * matrix elements are Mat3<double>
 * vector elements are Vec<3,double>
 */

#include <vector>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "Mat3.hpp"
#include "Vec.hpp"

template <typename Matrix, typename Vector, typename ResultVector>
ResultVector Matvec(const Matrix& A, const Vector& x)
{
  ResultVector r(x.size(),typename ResultVector::value_type(0));

  // get internal details from A
  auto& offsets = A.index1_data();
  auto& indices = A.index2_data();
  auto& values  = A.value_data();

  // loop over rows
  for (unsigned i=0; i<A.size1(); i++) {
    // loop over columns
    for (unsigned j=offsets[i]; j<offsets[i+1]; j++) {
      auto col = indices[j];
      auto temp = values[j]*x[col];
      r[i] = r[i]+temp;
    }
  }
  return r;
}

