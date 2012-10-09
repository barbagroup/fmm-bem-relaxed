#pragma once

/**
 * Storage class for all FMM options
 */

typedef enum {TOPDOWN, BOTTOMUP} TreeType;
typedef enum {FMM, TREECODE} EvaluatorType;

/** Class to define compile-time and run-time FMM options */
class FMM_options
{
public:
  bool symmetric;
  TreeType tree;
  EvaluatorType evaluator;
  double THETA;

  FMM_options() : symmetric(false), tree(TOPDOWN), evaluator(FMM), THETA(0.5) {}; 
};

