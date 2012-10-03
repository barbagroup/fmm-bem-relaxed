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

template <typename charge_type>
void chargesFromBodies(std::vector<charge_type>& charges, Bodies& bodies)
{
  for (size_t i=0; i<bodies.size(); i++)
  {
    charges[i] = charge_type(bodies[i].SRC);
  }
}

int main(int argc, char **argv)
{
  int numBodies = 100;
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

  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  Bodies bodies(numBodies);                                     // Define vector of bodies

  SphericalLaplaceKernel K(P);
  Dataset::cube(bodies,time(NULL));
  Bodies jbodies = bodies;                                               // Define vector of source bodies

  //fmm_plan plan = fmm_plan(K, bodies, opts);
  FMM_plan<SphericalLaplaceKernel> plan = FMM_plan<SphericalLaplaceKernel>(K, bodies, opts);

  std::vector<typename SphericalLaplaceKernel::charge_type> charges(numBodies);

  // get charges from initialised bodies
  chargesFromBodies(charges,bodies);

  //fmm_execute(plan, charges, jbodies);
  plan.execute(charges, jbodies);

  std::vector<typename SphericalLaplaceKernel::result_type> results(charges.size());
  plan.getResults(results.begin());

  // TODO: More elegant
  if (checkErrors) {
    const int numTargets = 100;                                  // Number of target points to be used for error eval
    bodies.resize(numTargets);
    Bodies test_bodies = bodies;

    for (B_iter B=test_bodies.begin(); B!=test_bodies.end(); ++B)
    {
      B->IBODY = B-test_bodies.begin();
      B->TRG = 0;
    }

    std::vector<typename SphericalLaplaceKernel::charge_type> source_charges(numBodies);
    std::vector<typename SphericalLaplaceKernel::point_type> source_points, test_points;
    std::vector<typename SphericalLaplaceKernel::result_type> test_results(numTargets);

    printf("Assembling test data\n");
    chargesFromBodies(source_charges, bodies);
    plan.bodies2points(jbodies, source_points);
    plan.bodies2points(test_bodies, test_points);
    printf("...done\n");

    printf("performing direct P2P\n");
    SimpleEvaluator<SphericalLaplaceKernel>::evalP2P(K,
                                                     source_points.begin(), source_points.end(), source_charges.begin(),
                                                     test_points.begin(), test_points.end(), test_results.begin());
    printf("...done\n");

    // TODO: Use a Direct class to make this more intuitive and accessible

    printf("checking errors\n");
    real diff1=0, diff2=0, norm1=0, norm2=0;
    auto fmm_it = results.begin();
    for (auto it=test_results.begin(); it!=test_results.end(); ++it, ++fmm_it)
    {
      printf("checking point..\n");
      std::cout << "Direct: " << *it << std::endl;
      std::cout << "FMM   : " << *fmm_it << std::endl;
      //printf("%lg\n",(*it)[0]);
      //printf("%lg\n",(*fmm_it)[0]);
      diff1 += ((*fmm_it)[0] - (*it)[0]) * ((*fmm_it)[0] - (*it)[0]);
      norm1 += (*it)[0] * (*it)[0];

      diff2 += ((*fmm_it)[1] - (*it)[1]) * ((*fmm_it)[1] - (*it)[1]);
      diff2 += ((*fmm_it)[2] - (*it)[2]) * ((*fmm_it)[2] - (*it)[2]);
      diff2 += ((*fmm_it)[3] - (*it)[3]) * ((*fmm_it)[3] - (*it)[3]);
      norm2 += (*it)[1]*(*it)[1];
      norm2 += (*it)[2]*(*it)[2];
      norm2 += (*it)[3]*(*it)[3];
    }
    printf("Error (pot) : %.4e\n",sqrt(diff1/norm1));
    printf("Error (acc) : %.4e\n",sqrt(diff2/norm2));
    

    #if 0
    B_iter B2 = bodies.begin();

    real diff1=0, diff2=0, norm1=0, norm2=0;
    for (B_iter B=test_bodies.begin(); B!=test_bodies.end(); ++B, ++B2)
    {
      diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);
      norm1 += B2->TRG[0] * B2->TRG[0];

      diff2 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);
      diff2 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);
      diff2 += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);
      norm2 += B2->TRG[1] * B2->TRG[1];
      norm2 += B2->TRG[2] * B2->TRG[2];
      norm2 += B2->TRG[3] * B2->TRG[3];
    }

    printf("Error (pot) : %.4e\n",sqrt(diff1/norm1));
    printf("Error (acc) : %.4e\n",sqrt(diff2/norm2));
    #endif
  }
}
