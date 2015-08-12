#pragma once

#include "FMM_plan.hpp"
#include "Preconditioner.hpp"
#include "GMRES.hpp"

FMMOptions& local_options()
{
  FMMOptions* opts = new FMMOptions;
  opts->local_evaluation = true;
  opts->lazy_evaluation = false;
  opts->set_mac_theta(0.5);
  opts->sparse_local = true;

  return *opts;
}


namespace Preconditioners {

/**
 * Preconditioner solving inner near-field problem
 */

template <class Plan, class Preconditioner>
class LocalInnerSolver
{
  typedef typename Plan::charge_type input_type;
  typedef typename Plan::result_type output_type;
  typedef std::vector<output_type> vector_type;
 private:
  Plan plan;
  vector_type& rhs;
  SolverOptions options;
  Preconditioner preconditioner;
  GMRESContext<output_type> context;

 public:
  //! Call inner solver from preconditioner interface
  template <class VecType>
  void operator()(VecType x, VecType& y) {
    // y = x;
    // printf("Calling Preconditioners::LocalInnerSolver\n");
    std::fill(y.begin(),y.end(),typename VecType::value_type(0.));
    GMRES(plan, y, x, options); // , preconditioner);
  }

  // construct from sources & targets
  template <typename Kernel, typename SourceVector, typename ResultVector>
  LocalInnerSolver(Kernel k, SourceVector& sources, ResultVector& RHS)
    : plan(k,sources,local_options()), rhs(RHS), preconditioner(k,sources.begin(),sources.end()), context(RHS.size(), 50) {
    //: plan(k,sources,local_options()), rhs(RHS), context(RHS.size(), 50) {

    options.residual= 1e-1;
    options.variable_p = false;
    options.max_iters = 1;

    context.output = false;
  }

  /*

  // need plan, RHS, options & preconditioner
  template <typename Vector>
  LocalInnerSolver(Plan& p, Vector& RHS, SolverOptions& opts, Preconditioner& P)
    : plan(p), rhs(RHS), options(opts), preconditioner(P), context(RHS.size(), opts.restart) {};

  // constructor using default preconditioner
  template <typename Vector>
  LocalInnerSolver(Plan& p, Vector& RHS, SolverOptions& opts)
    : LocalInnerSolver(p, RHS, opts, Preconditioners::Identity()) {};
  */
};

}; // end namespace Preconditioners
