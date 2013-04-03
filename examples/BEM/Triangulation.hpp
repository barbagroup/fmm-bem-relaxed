#pragma once

/**
 * holder for triangulation routines
 */

#include <cmath>
#include "Vec.hpp"
#include <fstream>

template <typename vector>
void normalise(vector& v)
{
  auto lens = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  v[0] /= lens;
  v[1] /= lens;
  v[2] /= lens;
}

struct triangle
{
  typedef struct triangle triangle_type;
  typedef Vec<3,double> vertex_type;
  vertex_type v0_;
  vertex_type v1_;
  vertex_type v2_;

  triangle() : v0_(0), v1_(0), v2_(0) {};
  triangle(vertex_type v0, vertex_type v1, vertex_type v2) : v0_(v0), v1_(v1), v2_(v2) {};
  // copy constructor
  triangle(const triangle& t) : v0_(t.v0_), v1_(t.v1_), v2_(t.v2_) {};
  // assignment
  triangle& operator=(const triangle& t)
  {
    v0_ = t.v0_;
    v1_ = t.v1_;
    v2_ = t.v2_;
    return *this;
  }
  /** subdivides into 4 triangles */
  std::vector<triangle_type> split() {
    std::vector<triangle_type> r(4);

    auto a = (v0_+v2_)*0.5;
    auto b = (v0_+v1_)*0.5;
    auto c = (v1_+v2_)*0.5;

    normalise(a);
    normalise(b);
    normalise(c);

    // assemble the 4 new triangles
    r[0] = triangle(v0_,b,a);
    r[1] = triangle(b,v1_,c);
    r[2] = triangle(a,b,c);
    r[3] = triangle(a,c,v2_);

    return r;
  }
};

std::vector<double> octahedron_vertices    = {1.,0.,0.,-1.,0.,0.,0.,1.,0.,0.,-1.,0.,0.,0.,1.,0.,0.,-1.};
std::vector<unsigned> octahedron_triangles = {0,4,2,2,4,1,1,4,3,3,4,0,0,2,5,2,1,5,1,3,5,3,0,5};

template <typename vertex_type>
void init_octahedron_vertices(std::vector<vertex_type>& vertices)
{
  for (unsigned i=0; i<6; i++) {
    vertices[i] = vertex_type(octahedron_vertices[i*3 + 0],
                              octahedron_vertices[i*3 + 1],
                              octahedron_vertices[i*3 + 2]);
  }
}

void init_octahedron_triangles(std::vector<triangle>& triangles)
{
  // get the initial vertices
  std::vector<triangle::vertex_type> initial_vertices(6);
  init_octahedron_vertices(initial_vertices);

  // construct the triangles
  for (unsigned i=0; i<8; i++) {
    triangles[i] = triangle(initial_vertices[octahedron_triangles[i*3 + 0]],
                            initial_vertices[octahedron_triangles[i*3 + 1]],
                            initial_vertices[octahedron_triangles[i*3 + 2]]);
  }
}

void divide_all(std::vector<triangle>& triangles)
{
  // take a copy of the triangles array & size
  auto triangles_copy = triangles;
  auto num_initial_triangles = triangles_copy.size();

  // resize output
  triangles.resize(num_initial_triangles*4);

  // loop over all triangles, splitting
  unsigned pos = 0;
  for (auto it=triangles_copy.begin(); it!=triangles_copy.end(); ++it, pos+=4) {
    auto split = it->split();

    for (unsigned i=0; i<4; i++) {
      triangles[pos + i] = split[i];
    }
  }
}

template <typename PanelType>
void create_unit_sphere(std::vector<PanelType>& panels, unsigned recursions = 2)
{
  std::vector<triangle> triangles(8);
  init_octahedron_triangles(triangles);

  for (unsigned i=0; i<recursions-1; i++) {
    divide_all(triangles);
  }

  printf("initialised %d triangles\n",(int)triangles.size());
  // now triangles contains the full triangulation -- convert to panel
  panels.resize(triangles.size());

  auto pit = panels.begin();
  for (auto it=triangles.begin(); it!=triangles.end(); ++it, ++pit) {
    *pit = PanelType(it->v0_,it->v1_,it->v2_);
  }

  // save the triangulation to test.vert, test.face
  std::ofstream face("test.face");
  std::ofstream vert("test.vert");

  int vnum = 1;
  for (auto it=triangles.begin(); it!=triangles.end();++it) {
    vert << it->v0_ << std::endl;
    vert << it->v1_ << std::endl;
    vert << it->v2_ << std::endl;
    face << vnum << ' ' << vnum+1 << ' ' << vnum+2 << std::endl;
    vnum+=3;
  }
}

