#pragma once

/**
 * Storage class for all FMM options
 */



/** Class to define compile-time and run-time FMM options */
class FMMOptions
{
public:
  typedef enum {TOPDOWN, BOTTOMUP} TreeType;
  typedef enum {FMM, TREECODE} EvaluatorType;

  bool symmetric;
  TreeType tree;
  EvaluatorType evaluator;
  double THETA;

  FMMOptions()
      : symmetric(false), tree(TOPDOWN), evaluator(FMM), THETA(0.5) {
  };
};

