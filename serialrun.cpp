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

template <typename Box>
void print_box(const Box& b, std::string padding = std::string()) {
  std::cout << padding << "Box " << b.index()
            << " (Level " << b.level() << ", Parent " << b.parent().index() << "): "
            << b.morton_index() << "    " << b.center() << "\n";

  padding.append(2,' ');
  if (!b.is_leaf()) {
    for (auto ci = b.child_begin(); ci != b.child_end(); ++ci)
      print_box(*ci, padding);
  } else {
    //for (auto ci = b.body_begin(); ci != b.body_end(); ++ci)
    //  std::cout << padding << "Point " << ci->index() << ": " << ci->morton_index() << "\t" << ci->point() << "\n";
  }
}

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}

int main(int argc, char **argv)
{
  int numBodies = 100;
  int P = 5;
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  bool checkErrors = true;
  FMM_options opts;

  // parse command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      i++;
      numBodies = atoi(argv[i]);
    } else if (strcmp(argv[i],"-P") == 0) {
      i++;
      P = atoi(argv[i]);
    } else if (strcmp(argv[i],"-theta") == 0) {
      i++;
      THETA = (double)atof(argv[i]);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    } else if (strcmp(argv[i],"-bottomup") == 0) {
      opts.tree = BOTTOMUP;
    } else {
      printf("[W]: Unknown command line arg: \"%s\"\n",argv[i]);
    }
  }

  SphericalLaplaceKernel K(P);
  typedef SphericalLaplaceKernel::point_type point_type;
  typedef SphericalLaplaceKernel::charge_type charge_type;
  typedef SphericalLaplaceKernel::result_type result_type;

  std::vector<point_type> points;
  for (int k = 0; k < numBodies; ++k)
    points.push_back(point_type(drand(), drand(), drand()));
  std::vector<point_type> jpoints = points;

  //fmm_plan plan = fmm_plan(K, bodies, opts);
  FMM_plan<SphericalLaplaceKernel> plan = FMM_plan<SphericalLaplaceKernel>(K, points, opts);
  print_box(plan.otree.root());

  std::vector<charge_type> charges(numBodies, charge_type(1));
  std::vector<charge_type> charges_copy = charges;

  //fmm_execute(plan, charges, jbodies);
  std::vector<result_type> result = plan.execute(charges, jpoints);

  // TODO: More elegant
  if (checkErrors) {
    std::vector<result_type> exact_result(numBodies, result_type(0));

    // TODO: Use a Direct class to make this more intuitive and accessible
    //SimpleEvaluator<SphericalLaplaceKernel>::evalP2P(K, test_bodies, jbodies);
    K.P2P(points.begin(), points.end(), charges_copy.begin(),
          jpoints.begin(), jpoints.end(),
          exact_result.begin());

    double diff1=0, diff2=0, norm1=0, norm2=0;
    int i=0;
    for (auto r1i = exact_result.begin(), r2i = result.begin();
             r1i != exact_result.end(); ++r1i, ++r2i) {

      result_type exact = *r1i;
      result_type r = *r2i;
      printf("[%03d] exact: %lg, FMM: %lg\n",i++,exact[0],r[0]);

      diff1 += (r[0] - exact[0]) * (r[0] - exact[0]);
      norm1 += exact[0] * exact[0];

      diff2 += (r[1] - exact[1]) * (r[1] - exact[1]);
      diff2 += (r[2] - exact[2]) * (r[2] - exact[2]);
      diff2 += (r[3] - exact[3]) * (r[3] - exact[3]);
      norm2 += exact[1] * exact[1];
      norm2 += exact[2] * exact[2];
      norm2 += exact[3] * exact[3];
    }


#if 0
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
#endif

    printf("Error (pot) : %.4e\n",sqrt(diff1/norm1));
    printf("Error (acc) : %.4e\n",sqrt(diff2/norm2));
  }
}
