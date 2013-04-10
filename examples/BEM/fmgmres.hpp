#pragma once

#include "Preconditioner.hpp"
#include "GMRES.hpp"

namespace Preconditioners {

template <class FMM_plan, class Preconditioner>
class FMGMRES
{
  // typedefs
  typedef typename FMM_plan::charge_type input_type;
  typedef typename FMM_plan::result_type output_type;
  typedef std::vector<output_type> vector_type;

 private:
  FMM_plan& plan;
  vector_type& rhs;
  const SolverOptions& options;
  Preconditioner& preconditioner;
  GMRESContext<output_type> context;

  /* interface to upload / download input / output vectors
  template <class VecType>
  void vector_upload(const VecType& in, vector_type& out) {
    auto out_it = out.begin();
    for (auto& it : in) {
      out_it++ = static_cast<output_type>(it);
    }
  }
  template <class VecType>
  void vector_download(const vector_type& in, VecType& out) {
    auto out_it = out.begin();
    for (auto& it : in) {
      out_it++ = static_cast<VecType::type>(it);
    }
  }
  */

 public:
  //! Call inner iteration from preconditioner interface
  // input / output vectors not necessarily same precision as input (worry about this later)
  // different precisions will require a conversion step (easy)
  template <class VecType>
  void operator()(VecType x, VecType& y) {
    GMRES(plan, y, x, options, preconditioner, context);
  }

  //! Need plan, RHS vector, Options & preconditioner (optional)
  template <typename Vector>
  FMGMRES(FMM_plan& Plan, Vector& RHS, SolverOptions& opts, Preconditioner& P)
    : plan(Plan), rhs(RHS), options(opts), preconditioner(P), context(RHS.size(),opts.restart) {};
  //! Constructor using default Preconditioner (by delagation)
  template <typename Vector>
  FMGMRES(FMM_plan& Plan, Vector& RHS, SolverOptions& opts)
    : FMGMRES(Plan, RHS, opts, Preconditioners::Identity()) {};

};

}; // end namespace preconditioners
