#pragma once

#include "EvaluatorBase.hpp"

#include "INITM.hpp"
#include "INITL.hpp"

#include "P2M.hpp"
#include "M2M.hpp"

template <typename Context, bool INITIALIZE_LOCAL>
struct EvalUpward : public EvaluatorBase<Context>
{
	void execute(Context& bc) const {
		auto& tree = bc.source_tree();

		// For the lowest level up to the highest level
		for (unsigned l = tree.levels()-1; l != 0; --l) {
			// For all boxes at this level
			auto b_end = tree.box_end(l);
			for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
				auto box = *bit;

				// TODO: initialize on-demand?
				INITM::eval(bc.kernel(), bc, box);
				if (INITIALIZE_LOCAL) INITL::eval(bc.kernel(), bc, box);

				if (box.is_leaf()) {
					// If leaf, make P2M calls
					P2M::eval(bc.kernel(), bc, box);
				} else {
					// If not leaf, then for all the children M2M
					auto c_end = box.child_end();
					for (auto cit = box.child_begin(); cit != c_end; ++cit)
						M2M::eval(bc.kernel(), bc, *cit, box);
				}
			}
		}
	}
};

template <typename Context, typename Options>
EvaluatorBase<Context>* make_upward(Context&, Options& opts) {
	if (opts.evaluator == FMMOptions::TREECODE)
		return new EvalUpward<Context, false>();
	if (opts.evaluator == FMMOptions::FMM)
		return new EvalUpward<Context, true>();
  return nullptr;
}
