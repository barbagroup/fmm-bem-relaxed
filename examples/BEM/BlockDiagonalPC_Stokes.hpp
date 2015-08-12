#pragma once

#include "FMM_plan.hpp"
#include "Preconditioner.hpp"
#include "GMRES_Stokes.hpp"



namespace Preconditioners {

/**
 * Preconditioner solving diagonal inner near-field problem
 */

template <class Plan>
class BlockDiagonal
{
  typedef typename Plan::charge_type input_type;
  typedef typename Plan::result_type output_type;
  typedef std::vector<output_type> vector_type;
 private:
  Plan plan;
  SolverOptions options;
  GMRESContext<output_type> context;
  Preconditioners::Identity M;

 public:
  //! Call inner solver from preconditioner interface
  template <class VecType>
  void operator()(VecType x, VecType& y) {
    // y = x;
    // printf("Calling Preconditioners::LocalInnerSolver\n");
    std::fill(y.begin(),y.end(),typename VecType::value_type(0.));
    GMRES(plan, y, x, options, M, context);
    context.reset();
  }

  // construct from sources & targets
  template <typename Kernel, typename SourceVector>
  BlockDiagonal(Kernel k, SourceVector& sources)
    : plan(k,sources,BlockDiagonal::local_options()), context(sources.size(), 50), M() {

    options.residual= 1e-1;
    options.variable_p = false;
    options.max_iters = 1;

    context.output = false;
  }

  static FMMOptions& local_options()
  {
    FMMOptions* opts = new FMMOptions;
    opts->local_evaluation = false;
    opts->lazy_evaluation = false;
    opts->set_mac_theta(0.5);
    opts->sparse_local = true;
    opts->block_diagonal = true;
    opts->set_max_per_box(100);

    return *opts;
  }

};

}; // end namespace Preconditioners
