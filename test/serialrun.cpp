/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <FMM_plan.hpp>
#include <Dataset.hpp>
#include <SphericalLaplaceKernel.hpp>

int main(int argc, char **argv)
{
  int numBodies = 10000;
  int P = 8;
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  bool checkErrors = true;
  FMM_options opts;

  // parse command line args
  for (int i=1; i<argc; i++)
  {
    if (strcmp(argv[i],"-N")==0)
    {
      i++;
      numBodies = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-P")==0)
    {
      i++;
      P = atoi(argv[i]);
    }
    else if (strcmp(argv[i],"-theta")==0)
    {
      i++;
      THETA = (double)atof(argv[i]);
    }
    else if (strcmp(argv[i],"-nocheck")==0)
    {
      checkErrors = false;
    }
    else if (strcmp(argv[i],"-bottomup")==0)
    {
      opts.tree = BOTTOMUP;
    }
    else
    {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
    }
  }

  const int numTarget = 100;                                    // Number of target points to be used for error eval
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  Bodies bodies(numBodies);                                     // Define vector of bodies

  SphericalLaplaceKernel K(P);
  Dataset::cube(bodies,time(NULL));
  Bodies jbodies = bodies;                                               // Define vector of source bodies

  FMM_plan<SphericalLaplaceKernel> FMM(K,bodies,opts);
  FMM.execute(jbodies);
  if (checkErrors)
    FMM.checkError(bodies,jbodies);
}
