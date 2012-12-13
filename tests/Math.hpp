#pragma once

/** Compute the L1 error of v1 and v2
 * E = 1/N * sum_k |v1[k] - v2[k]|
 *
 * @pre The class Vec needs .size() and .op[]
 */
template <typename Vec>
double l1_error(const Vec& v1, const Vec& v2) {
  assert(v1.size() == v2.size());
  unsigned N = v1.size();
  double e = 0;
  for (unsigned k = 0; k < N; ++k)
    e += fabs(v1[k] - v2[k]);
  return e / N;
}


/** Compute the L1 relative error of v with respect to ve
 * E = 1/N * sum_k |ve[k] - v[k]| / |ve[k]|
 *
 * @pre The class Vec needs .size() and .op[]
 */
template <typename Vec>
double l1_rel_error(const Vec& v, const Vec& ve) {
  assert(v.size() == ve.size());
  unsigned N = v.size();
  double e = 0;
  double ev = 0;
  for (unsigned k = 0; k < N; ++k) {
    e += fabs(ve[k] - v[k]);
    ev += fabs(ve[k]);
  }
  return e / (N*ev);
}

/** Compute the L2 error of v1 and v2
 * E = sqrt( sum_k (ve[k] - v[k])^2 )
 *
 * @pre The class Vec needs .size() and .op[]
 */
template <typename Vec>
double l2_error(const Vec& v1, const Vec& v2) {
  assert(v1.size() == v2.size());
  unsigned N = v1.size();
  double e2 = 0;
  for (unsigned k = 0; k < N; ++k) {
    double e = v1[k] - v2[k];
    e2 += e * e;
  }
  return sqrt(e2);
}


/** Compute the L2 relative error of v with respect to ve
 * E = sqrt( sum_k (ve[k] - v[k])^2 / ve[k]^2 )
 *
 * @pre The class Vec needs .size() and .op[]
 */
template <typename Vec>
double l2_rel_error(const Vec& v, const Vec& ve) {
  assert(v.size() == ve.size());
  unsigned N = v.size();
  double e2 = 0;
  double ve2 = 0;
  for (unsigned k = 0; k < N; ++k) {
    double e = ve[k] - v[k];
    e2 += e * e;
    ve2 += ve[k] * ve[k];
  }
  return sqrt(e2 / ve2);
}

/** Compute the max error of v1 and v2
 * E = max_k |ve[k] - v[k]|
 *
 * @pre The class Vec needs .size() and .op[]
 */
template <typename Vec>
double max_error(const Vec& v1, const Vec& v2) {
  assert(v1.size() == v2.size());
  unsigned N = v1.size();
  double e = 0;
  for (unsigned k = 0; k < N; ++k)
    e = max(e, fabs(v1[k] - v2[k]));
  return e;
}


/** Compute the max relative error of v with respect to ve
 * E = max_k |ve[k] - v[k]| / |ve[k]|
 *
 * @pre The class Vec needs .size() and .op[]
 */
template <typename Vec>
double max_rel_error(const Vec& v, const Vec& ve) {
  assert(v.size() == ve.size());
  unsigned N = v.size();
  double e = 0;
  for (unsigned k = 0; k < N; ++k)
    e = max(e, fabs(ve[k] - v[k]) / fabs(ve[k]));
  return e;
}
