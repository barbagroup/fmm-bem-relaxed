#pragma once

/**
 * Storage class for all FMM options
 */


/** Class to define compile-time and run-time FMM options */
class FMMOptions
{
public:
	bool lazy_evaluation; // use lazy evaluation of multipole / local expansions / translations?
  bool local_evaluation; // use only local evaluation (for preconditioners)

	//! Evaluation type
	enum EvalType {FMM, TREECODE};
	EvalType evaluator;


	struct DefaultMAC {
		double theta_;
		DefaultMAC(double theta) : theta_(theta) {}

		template <typename BOX>
		bool operator()(const BOX& b1, const BOX& b2) const {
			double r0_norm = norm(b1.center() - b2.center());
			return r0_norm * theta_ > b1.radius() + b2.radius();
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
		  evaluator(FMM),
		  MAC_(DefaultMAC(0.5)),
		  NCRIT_(126),
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Get the FMMOptions from command line arguments
 */
FMMOptions get_options(int argc, char** argv) {
	FMMOptions opts;
	opts.set_mac_theta(0.5);    // Default multipole acceptance criteria
	opts.set_max_per_box(64);   // Default NCRIT

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
