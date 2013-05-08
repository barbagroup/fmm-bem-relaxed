#pragma once

/*
 * Main class to hold shared things for (semi-) analytical integration routines
 */

namespace AnalyticalIntegral
{
  // type of equation
  typedef enum { LAPLACE, YUKAWA, HELMHOLTZ, STOKES } equation;

}; // end namespace AnalyticalIntegral
