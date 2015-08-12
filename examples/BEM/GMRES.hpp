#pragma once

/**
 * New implementation of GMRES / FGMRES based on abstractions to Matvec and storage
 */

#include <cstdio>
#include <cmath>
#include <vector>
#include <type_traits>

#include "Preconditioner.hpp"
#include "SolverOptions.hpp"
#include "BLAS.hpp"

/** default case */
template <typename T, typename ...Args>
struct ChargeSize {
  static int size() { return 1; }
};

template <typename T>
struct ChargeSize<T> {
  static int size() { return 1; }
};
template <int N, typename T2>
struct ChargeSize<Vec<N,T2>> {
  static int size() { return N; }
};

template <>
struct ChargeSize<Vec<3,double>> {
  static int size() { return 3;}
};

//! Context to store all GMRES temporary vectors -- shared between calls
template <typename T>
struct GMRESContext
{
 public:
  // vector quantities
  std::vector<T> w, V0, z, s, cs, sn;

  // matrix quantities
  Matrix<T> V, H;

  bool output;

  //! constructor
  GMRESContext(const unsigned N=0, const unsigned R=50) : output(true)
//     : w(N), V0(N), z(N), s(R+1), cs(R), sn(R), V(N,R+1), H(R+1,R) {};
  {
    unsigned factor = ChargeSize<T>::size();
    w  = std::vector<T>(factor*N);
    V0 = std::vector<T>(factor*N);
    z  = std::vector<T>(factor*N);
    s  = std::vector<T>(R+1);
    cs = std::vector<T>(R);
    sn = std::vector<T>(R);
    V  = Matrix<T>(factor*N,R+1);
    H  = Matrix<T>(R+1,R);
  }
};

//! Context for FGMRES temporary vectors
template <typename T>
struct FGMRESContext : public GMRESContext<T>
{
 public:
  Matrix<T> Z;

  //! constructor
  FGMRESContext(const unsigned N, const unsigned R)
    : GMRESContext<T>(N,R), Z(N,R+1) {};
};

/** GMRES implementation
 * requires Matvec object with execute(std::vector<T>&, unsigned) signature
 * Matvec must also define charge and result types
 */

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

//! pass with neither preconditioner or context
template <typename Matvec>
void GMRES(Matvec& MV,
           std::vector<typename Matvec::charge_type>& x,
           std::vector<typename Matvec::result_type>& b,
           const SolverOptions& opts)
{
  typedef typename Matvec::result_type result_type;
  GMRESContext<result_type> context(x.size(), opts.restart);
  GMRES(MV,x,b,opts,Preconditioners::Identity(),context);
}

//! pass with preconditioner but no context
template <typename Matvec, typename Preconditioner>
void GMRES(Matvec& MV,
           std::vector<typename Matvec::charge_type>& x,
           std::vector<typename Matvec::result_type>& b,
           const SolverOptions& opts, const Preconditioner& M)
{
  typedef typename Matvec::result_type result_type;
  GMRESContext<result_type> context(x.size(), opts.restart);
  GMRES(MV,x,b,opts,M,context);
}

template <typename Matvec, typename Preconditioner, typename SolverContext>
void GMRES(Matvec& MV,
           std::vector<typename Matvec::charge_type>& x,
           std::vector<typename Matvec::result_type>& b,
           const SolverOptions& opts, const Preconditioner& M,
           SolverContext& context)
{
  typedef typename Matvec::charge_type charge_type;
  typedef typename Matvec::result_type result_type;

  // get sizes of charges and results
  //int charge_size = ChargeSize<charge_type>::size();
  //int result_size = ChargeSize<result_type>::size();

  const int R = opts.restart;
  const int max_iters = opts.max_iters;
  result_type beta = 0.;

  int i; // , j;
  int iter = 0;
  // allocate the workspace
  result_type resid = 0;

  // scale residual by ||b||
  auto normb = norm(b.begin(),b.end());

  // reference to the kernel used
  auto& K = MV.kernel();

  // outer (restart) loop
  do {
    // dot product of A*x -- FMM call
    std::fill(context.w.begin(),context.w.end(),0);
    context.w = MV.execute(x); // V(0) = A*x
    // V(0) = V(0) - b
    blas::axpy(b,context.w,-1.);
    beta = blas::nrm2(context.w); // beta = ||V(0)||
    blas::scal(context.w,-1./beta); // V(0) = -V(0)/beta
    // w at this point contains V(0)
    // V(:,0) = V(0) { V(:,0) = w }
    context.V.set_column(0,context.w);
    // set S = 0

    context.s[0] = beta;
    i = -1;
    resid = context.s[0]/normb;

    // inner loop
    do {
      ++i;
      ++iter;

      // set p for this iteration
      int p = std::max(1u,opts.predict_p(fabs(resid)));
      K.set_p(p);

      // perform w = A*x
      std::fill(context.V0.begin(),context.V0.end(),0.);
      M(context.V.column(i),context.z);
      context.w = MV.execute(context.z);

      for (int k=0; k<=i; k++) {
        auto V_col = context.V.column(k);
        context.H(k,i) = blas::dotc(context.w,V_col);
        // V(i+1) -= H(k,i)*V(k)
        blas::axpy(V_col,context.w,-context.H(k,i));
      }

      // H(i+1,i) = nrm2(w)
      context.H(i+1,i) = blas::nrm2(context.w);
      // V(i+1) = V(i+1) / H(i+1,i)
      auto w_temp = context.w;
      blas::scal(w_temp,1./context.H(i+1,i));
      // V(:,i+1) = w
      context.V.set_column(i+1,w_temp);

      PlaneRotation(context.H,context.cs,context.sn,context.s,i);

      resid = context.s[i+1]/normb;

      if (fabs(resid) < opts.residual) break;

      if (context.output)
        printf("it: %03d, res: %.3e, fmm_req_p: %01d\n",iter,fabs(resid),p);

    } while(i+1 < R && i+1 <= max_iters  && fabs(resid) > opts.residual);

    // solve upper triangular system in place
    for (int j=i; j >= 0; j--) {
      context.s[j] /= context.H(j,j);
      // s[0:j] = s[0:j] - s[j]*H(0:j,j)
      for (int k=j-1; k>=0; k--) {
        context.s[k] -= context.H(k,j)*context.s[j];
      }
    }

    // Update the solution
    for (int j=0; j<=i; j++) {
      // x =x + x[j]*V(:,j)
      M(context.V.column(j),context.z);
      blas::axpy(context.z,x,context.s[j]);
    }
    if (context.output)
      if (iter % 10 == 0) printf("it: %04d, residual: %.3e\n",iter,(double)fabs(resid));

  } while (fabs(resid) > opts.residual && iter < opts.max_iters);
  // } while (iter < opts.max_iters);

  if (context.output)
    printf("Final residual: %.4e, after %d iterations\n",fabs(resid),iter);
}

template <typename Matvec>
void FGMRES(Matvec& MV,
            std::vector<typename Matvec::charge_type>& x,
            std::vector<typename Matvec::result_type>& b,
            const SolverOptions& opts)
{
  Preconditioners::Identity M;
  FGMRES(MV,x,b,opts,M);
}

//! pass with preconditioner but no context
template <typename Matvec, typename Preconditioner>
void FGMRES(Matvec& MV,
            std::vector<typename Matvec::charge_type>& x,
            std::vector<typename Matvec::result_type>& b,
            const SolverOptions& opts, Preconditioner& M)
{
  typedef typename Matvec::result_type result_type;
  FGMRESContext<result_type> context(x.size(), opts.restart);
  FGMRES(MV,x,b,opts,M,context);
}

template <typename Matvec, typename Preconditioner, typename SolverContext>
void FGMRES(Matvec& MV,
            std::vector<typename Matvec::charge_type>& x,
            std::vector<typename Matvec::result_type>& b,
            const SolverOptions& opts, Preconditioner& M,
            SolverContext& context)
{
  typedef typename Matvec::charge_type charge_type;
  typedef typename Matvec::result_type result_type;

  const int R = opts.restart;
  const int max_iters = opts.max_iters;
  result_type beta = 0.;

  int i; // , j;
  int iter = 0;
  // allocate the workspace
  result_type resid = 0;

  // scale residual by ||b||
  auto normb = norm(b.begin(),b.end());

  auto& K = MV.kernel();

  // outer (restart) loop
  do {
    // dot product of A*x -- FMM call
    context.w = MV.execute(x); // V(0) = A*x
    // V(0) = V(0) - b
    blas::axpy(b,context.w,-1.);
    beta = blas::nrm2(context.w); // beta = ||V(0)||
    blas::scal(context.w,-1./beta); // V(0) = -V(0)/beta
    // w at this point contains V(0)
    // V(:,0) = V(0) { V(:,0) = w }
    context.V.set_column(0,context.w);
    // set S = 0

    context.s[0] = beta;
    i = -1;
    resid = context.s[0]/normb;

    // inner loop
    do {
      ++i;
      ++iter;

      // set p for this iteration
      int p = opts.predict_p(fabs(resid));
      K.set_p(p);

      // perform w = A*x
      std::fill(context.V0.begin(),context.V0.end(),0.);
      M(context.V.column(i),context.z);
      context.Z.set_column(i,context.z);
      context.w = MV.execute(context.z);

      for (int k=0; k<=i; k++) {
        auto V_col = context.V.column(k);
        context.H(k,i) = blas::dotc(context.w,V_col);
        // V(i+1) -= H(k,i)*V(k)
        blas::axpy(V_col,context.w,-context.H(k,i));
      }

      // H(i+1,i) = nrm2(w)
      context.H(i+1,i) = blas::nrm2(context.w);
      // V(i+1) = V(i+1) / H(i+1,i)
      auto w_temp = context.w;
      blas::scal(w_temp,1./context.H(i+1,i));
      // V(:,i+1) = w
      context.V.set_column(i+1,w_temp);

      PlaneRotation(context.H,context.cs,context.sn,context.s,i);

      resid = context.s[i+1]/normb;

      if (fabs(resid) < opts.residual) break;

      if (context.output)
        printf("it: %03d, res: %.3e, fmm_req_p: %01d\n",iter,fabs(resid),p);

    } while(i+1 < R && i+1 <= max_iters && fabs(resid) > opts.residual);

    // solve upper triangular system in place
    for (int j=i; j >= 0; j--) {
      context.s[j] /= context.H(j,j);
      // s[0:j] = s[0:j] - s[j]*H(0:j,j)
      for (int k=j-1; k>=0; k--) {
        context.s[k] -= context.H(k,j)*context.s[j];
      }
    }

    // Update the solution
    for (int j=0; j<=i; j++) {
      // x =x + x[j]*V(:,j)
      blas::axpy(context.Z.column(j),x,context.s[j]);
    }

    if (context.output)
      if (iter % 10 == 0) printf("it: %04d, residual: %.3e\n",iter,(double)fabs(resid));

  } while (fabs(resid) > opts.residual && iter < opts.max_iters);

  if (context.output)
    printf("Final residual: %.4e, after %d iterations\n",fabs(resid),iter);
}
