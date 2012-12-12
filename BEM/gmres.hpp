#pragma once

/**
 * gmres.hpp
 *
 * implementation of GMRES krylov solver
 * Code modified from Cusp's cusp::krylov::gmres
 */

#include <cstdio>
#include <vector>
#include <cmath>

/** Very basic matrix type */
template <typename T>
class Matrix
{
 private:
  unsigned rows_, cols_;
  // 1D storage, column-major
  std::vector<T> vals_;

 public:
  Matrix() : rows_(0), cols_(0), vals_(0) {};
  Matrix(int r, int c) : rows_(r), cols_(c) {
    vals_.resize(rows_*cols_);
  };
  Matrix(const Matrix& M) : rows_(M.rows_), cols_(M.cols_), vals_(M.vals_) {};

  unsigned rows() { return rows_; };
  unsigned cols() { return cols_; };

  std::vector<T> column(int c) {
    std::vector<T> r(rows_,0);
    for (unsigned i=0; i<rows_; i++) {
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

template <typename T>
void scal(std::vector<T>& vec, T s)
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

}; // end namespace blas

template <typename Iter>
typename Iter::value_type norm(Iter start, Iter end)
{
  typename Iter::value_type sum = 0;

  for ( ; start != end; ++start) sum += (*start) * (*start);

  return std::sqrt(sum);
}

template <typename T>
void ApplyPlaneRotation(T& dx, T& dy, T& cs, T& sn)
{
  T temp = cs * dx + sn *dy;
  dy = -sn*dx+cs*dy;
  dx = temp;
}

template <typename T>
void GeneratePlaneRotation(T& dx, T&dy, T& cs, T& sn)
{
  if(dy == T(0.0)){
    cs = 1.0;
    sn = 0.0;
  }else if (fabs(dy) > fabs(dx)) {
    T tmp = dx / dy;
    sn = T(1.0) / sqrt(T(1.0) + tmp*tmp);
    cs = tmp*sn;            
  }else {
    T tmp = dy / dx;
    cs = T(1.0) / sqrt(T(1.0) + tmp*tmp);
    sn = tmp*cs;
  }
}

template <class Matrix, typename T>
void PlaneRotation(Matrix& H, T& cs, T& sn, T& s, int i)
{
  for (int k = 0; k < i; k++){
    ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
  }
  GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
  ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
  ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
}

template <typename FMM, typename SolverOptions>
void fmm_gmres(FMM& fmm,
               std::vector<typename FMM::charge_type>& x,
               std::vector<typename FMM::result_type>& b,
               const SolverOptions& opts)
{
  const int N = x.size();
  const int R = opts.restart;
  const int max_iters = opts.max_iters;
  double beta = 0.;

  int i; // , j;
  int iter = 0;
  // allocate the workspace
  double resid = 0;
  std::vector<typename FMM::result_type> w(N);
  std::vector<double> V0(N);
  Matrix<double> V(N,R+1);
  Matrix<double> H(R+1,R);
  std::vector<double> s(R+1);
  std::vector<double> cs(R);
  std::vector<double> sn(R);
  // bool converged = false;

  // scale residual by ||b||
  auto normb = norm(b.begin(),b.end());

  // main loop
  do {
    // dot product of A*x -- FMM call
    w = fmm.execute(x); // V(0) = A*x
    // V(0) = V(0) - b
    blas::axpy(b,w,-1.);
    // V(0) = M*V(0) -- preconditioner step if desired
    beta = blas::nrm2(w); // beta = ||V(0)||
    blas::scal(w,-1./beta); // V(0) = -V(0)/beta
    // V(:,0) = V(0) { V(:,0) = w }
    V.set_column(0,w);
    // set S = 0

    s[0] = beta;
    i = -1;
    resid = s[0]/normb;


    // inner loop
    do {
      ++i;
      ++iter;

      // perform w = A*x
      std::fill(V0.begin(),V0.end(),0.);
      V0 = fmm.execute(V.column(i));
      w = V0; // instead of preconditioner step

      for (int k=0; k<=i; k++) {
        auto V_col = V.column(k);
        H(k,i) = blas::dotc(w,V_col);
        // V(i+1) -= H(k,i)*V(k)
        blas::axpy(V_col,w,-H(k,i));
      }

      // H(i+1,i) = nrm2(w)
      H(i+1,i) = blas::nrm2(w);
      // V(i+1) = V(i+1) / H(i+1,i)
      auto w_temp = w;
      blas::scal(w_temp,1./H(i+1,i));
      // V(:,i+1) = w
      V.set_column(i+1,w_temp);

      PlaneRotation(H,cs,sn,s,i);

      resid = s[i+1]/normb;

      if (fabs(resid) < opts.residual) break;
      printf("it: %03d, res: %.3e\n",iter,fabs(resid));

    } while(i+1 < R && i+1 <= max_iters && fabs(resid) > opts.residual);

    // solve upper triangular system in place
    for (int j=i; j >= 0; j--) {
      s[j] /= H(j,j);
      // s[0:j] = s[0:j] - s[j]*H(0:j,j)
      for (int k=j-1; k>=0; k--) {
        s[k] -= H(k,j)*s[j];
      }
    }

    // Update the solution
    for (int j=0; j<=i; j++) {
      // x =x + x[j]*V(:,j)
      blas::axpy(V.column(j),x,s[j]);
    }
    if (iter % 10 == 0) printf("it: %04d, residual: %.3e\n",iter,(double)fabs(resid));

  } while (fabs(resid) > opts.residual && iter < opts.max_iters);
  printf("Final residual: %.4e, after %d iterations\n",fabs(resid),iter);
}

template <typename kernel, typename source_type, typename charge_type, typename result_type, typename SolverOptions>
void direct_gmres(kernel& k, std::vector<source_type>& panels,
                  std::vector<charge_type>& x,
                  std::vector<result_type>& b,
                  const SolverOptions& opts)
{
  const int N = x.size();
  const int R = opts.restart;
  const int max_iters = opts.max_iters;
  double beta = 0.;

  int i; // , j;
  int iter = 0;
  // allocate the workspace
  double resid = 0;
  std::vector<result_type> w(N);
  std::vector<double> V0(N);
  Matrix<double> V(N,R+1);
  Matrix<double> H(R+1,R);
  std::vector<double> s(R+1);
  std::vector<double> cs(R);
  std::vector<double> sn(R);
  // bool converged = false;

  // main loop
  do {
    // dot product of A*x -- FMM call
    // w = fmm.execute(x); // V(0) = A*x
    Direct::matvec(k,panels,x,b);
    //printf("V(0) = A*x\n");
    // V(0) = V(0) - b
    blas::axpy(b,w,-1.);
    //printf("V(0) = -V(0) / beta\n");
    // V(0) = M*V(0) -- preconditioner step if desired
    beta = blas::nrm2(w); // beta = ||V(0)||
    blas::scal(w,-1./beta); // V(0) = -V(0)/beta
    // V(:,0) = V(0) { V(:,0) = w }
    V.set_column(0,w);
    // set S = 0

    s[0] = beta;
    std::cout << "beta " << beta << std::endl;
    i = -1;
    resid = s[0];


    // inner loop
    do {
      ++i;
      ++iter;

      // perform w = A*x
      // V0 = fmm.execute(w);
      std::fill(V0.begin(),V0.end(),0.);
      Direct::matvec(k,panels,V.column(i),V0);
      w = V0; // instead of preconditioner step


      for (int k=0; k<=i; k++) {
        // H(k,i) = dotc(w,V(:,k));
        auto V_col = V.column(k);
        H(k,i) = blas::dotc(w,V_col);
        // V(i+1) -= H(k,i)*V(k)
        blas::axpy(V_col,w,-H(k,i));
      }

      // H(i+1,i) = nrm2(w)
      H(i+1,i) = blas::nrm2(w);
      // V(i+1) = V(i+1) / H(i+1,i)
      auto w_temp = w;
      blas::scal(w_temp,1./H(i+1,i));
      // V(:,i+1) = w
      V.set_column(i+1,w_temp);

      PlaneRotation(H,cs,sn,s,i);

      resid = s[i+1];

      if (fabs(resid) < opts.residual) break;
      printf("it: %03d, res: %.3e\n",iter,fabs(resid));

    } while(i+1 < R && i+1 <= max_iters && fabs(resid) > opts.residual);

    // solve upper triangular system in place
    for (int j=i; j >= 0; j--) {
      s[j] /= H(j,j);
      // s[0:j] = s[0:j] - s[j]*H(0:j,j)
      for (int k=j-1; k>=0; k--) {
        s[k] -= H(k,j)*s[j];
      }
    }

    // Update the solution
    for (int j=0; j<=i; j++) {
      // x =x + x[j]*V(:,j)
      blas::axpy(V.column(j),x,s[j]);
    }
    if (iter % 10 == 0) printf("it: %04d, residual: %.3e\n",iter,(double)fabs(resid));

  } while (fabs(resid) > opts.residual && iter < opts.max_iters);
  printf("Final residual: %.4e, after %d iterations\n",fabs(resid),iter);
}

template <typename Mat, typename SolverOptions>
void gmres(Mat& A,
           std::vector<double>& x,
           std::vector<double>& b,
           const SolverOptions& opts)
{
  typedef double source_type;
  typedef double charge_type;
  typedef double result_type;

  const int N = x.size();
  const int R = opts.restart;
  const int max_iters = opts.max_iters;
  double beta = 0.;

  int i; // , j;
  int iter = 0;
  // allocate the workspace
  double resid = 0;
  std::vector<result_type> w(N);
  std::vector<double> V0(N);
  Matrix<double> V(N,R+1);
  Matrix<double> H(R+1,R);
  std::vector<double> s(R+1);
  std::vector<double> cs(R);
  std::vector<double> sn(R);
  // bool converged = false;

  // main loop
  do {
    iter++;
    // dot product of A*x -- FMM call
    // w = FMM.execute(x); // V(0) = A*x
    blas::matvec(A,x,w);
    //printf("V(0) = A*x\n");
    // V(0) = V(0) - b
    blas::axpy(b,w,-1.);
    //printf("V(0) = -V(0) / beta\n");
    // V(0) = M*V(0) -- preconditioner step if desired
    beta = blas::nrm2(w); // beta = ||V(0)||
    blas::scal(w,-1./beta); // V(0) = -V(0)/beta
    // V(:,0) = V(0) { V(:,0) = w }
    V.set_column(0,w);
    // set S = 0
    for (unsigned l=0; l<s.size(); l++) s[l] = 0.;

    s[0] = beta;
    i = -1;
    resid = s[0];

    // inner loop
    do {
      ++i;

      // perform w = A*x
      // V0 = FMM.execute(x);
      blas::matvec(A,w,V0);
      w = V0; // instead of preconditioner step

      for (int k=0; k<=i; k++) {
        // H(k,i) = dotc(w,V(:,k));
        auto V_col = V.column(k);
        H(k,i) = blas::dotc(w,V_col);
        // V(i+1) -= H(k,i)*V(k)
        blas::axpy(V_col,w,-H(k,i));
      }

      // H(i+1,i) = nrm2(w)
      H(i+1,i) = blas::nrm2(w);
      // V(i+1) = V(i+1) / H(i+1,i)
      blas::scal(w,1./H(i+1,i));
      // V(:,i+1) = w
      V.set_column(i+1,w);

      PlaneRotation(H,cs,sn,s,i);

      resid = s[i+1];

      if (fabs(resid) < opts.residual) break;

    } while(i+1 < R && i+1 <= max_iters && fabs(resid) > opts.residual);

    // solve upper triangular system in place
    for (int j=i; j >= 0; j--) {
      s[j] /= H(j,j);
      // s[0:j] = s[0:j] - s[j]*H(0:j,j)
      for (int k=j-1; k>=0; k--) {
        s[k] -= H(k,j)*s[j];
      }
    }

    // Update the solution
    for (int j=0; j<=i; j++) {
      // x =x + x[j]*V(:,j)
      blas::axpy(V.column(j),x,s[j]);
    }
    if (iter % 10 == 0) printf("it: %04d, residual: %.3e\n",iter,(double)fabs(resid));

  } while (fabs(resid) > opts.residual && iter < opts.max_iters);
  printf("Final residual: %.4e, after %d iterations\n",fabs(resid),iter);
}
