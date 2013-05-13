#pragma once
/**
 * Storage class for all FMM options
 */

#include "Vec.hpp"

/** Class to define compile-time and run-time FMM options */
class FMMOptions
{
public:
	bool lazy_evaluation; // use lazy evaluation of multipole / local expansions / translations?
  bool local_evaluation; // use only local evaluation (for preconditioners)
  bool sparse_local; // only local eval using sparse matrix
  bool block_diagonal; // only diagonal eval, using sparse matrix

	//! Evaluation type
	enum EvalType {FMM, TREECODE};
	EvalType evaluator;

	struct DefaultMAC {
		double theta_;
		DefaultMAC(double theta) : theta_(theta) {}

		template <typename BOX>
		bool operator()(const BOX& b1, const BOX& b2) const {
			double r0_normSq = normSq(b1.center() - b2.center());
      double rhs = (b1.radius() + b2.radius()) / theta_;
			return r0_normSq > rhs*rhs;
		}
	};

	// TODO: Generalize type?
	DefaultMAC MAC_;
	unsigned NCRIT_;

	// DEBUGGING FLAGS
	bool printTree;

	FMMOptions()
		: lazy_evaluation(false),
      local_evaluation(false),
      sparse_local(false),
      block_diagonal(false),
		  evaluator(FMM),
		  MAC_(DefaultMAC(0.5)),
		  NCRIT_(64),
		  printTree(false) {
	};

	void set_mac_theta(double theta) {
		MAC_ = DefaultMAC(theta);
	}

	DefaultMAC MAC() {
		return MAC_;
	}

	void set_max_per_box(unsigned ncrit) {
		NCRIT_ = ncrit;
	}

	unsigned max_per_box() const {
		return NCRIT_;
	}

	void print_tree(bool v) { printTree = v; }
	bool print_tree() const { return printTree; }
};


#include <cstdio>
#include <cstdlib>
#include <cstring>

/** Get the FMMOptions from command line arguments
 */
FMMOptions get_options(int argc, char** argv) {
	FMMOptions opts = FMMOptions();

	// parse command line args
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i],"-theta") == 0) {
			i++;
			opts.set_mac_theta((double)atof(argv[i]));
		} else if (strcmp(argv[i],"-eval") == 0) {
			i++;
			if (strcmp(argv[i],"FMM") == 0) {
				opts.evaluator = FMMOptions::FMM;
			} else if (strcmp(argv[i],"TREE") == 0) {
				opts.evaluator = FMMOptions::TREECODE;
			} else {
				printf("[W]: Unknown evaluator type: \"%s\"\n",argv[i]);
			}
		} else if (strcmp(argv[i],"-lazy_eval") == 0) {
			opts.lazy_evaluation = true;
		} else if (strcmp(argv[i],"-ncrit") == 0) {
			i++;
			opts.set_max_per_box((unsigned)atoi(argv[i]));
		} else if (strcmp(argv[i],"-printtree") == 0) {
			opts.print_tree(true);
		}
	}

	return opts;
}
