#pragma once

/**
 * holder for triangulation routines
 */

#include <cmath>
#include "Vec.hpp"
#include <fstream>
#include <cstdlib>

namespace Triangulation {

struct triangle
{
  typedef struct triangle triangle_type;
  typedef Vec<3,double> vertex_type;
  vertex_type v0_;
  vertex_type v1_;
  vertex_type v2_;

  triangle() {};
  triangle(vertex_type v0, vertex_type v1, vertex_type v2)
      : v0_(v0), v1_(v1), v2_(v2) {}
  // copy constructor
  triangle(const triangle& t)
      : v0_(t.v0_), v1_(t.v1_), v2_(t.v2_) {}
  // assignment
  triangle& operator=(const triangle& t) {
    v0_ = t.v0_;
    v1_ = t.v1_;
    v2_ = t.v2_;
    return *this;
  }
  /** subdivides into 4 triangles */
  std::vector<triangle_type> split() {
    std::vector<triangle_type> r(4);

    vertex_type a = (v0_+v2_)*0.5;
    vertex_type b = (v0_+v1_)*0.5;
    vertex_type c = (v1_+v2_)*0.5;

    a /= norm(a);
    b /= norm(b);
    c /= norm(c);

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
void UnitSphere(std::vector<PanelType>& panels, unsigned recursions = 2)
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
    vert << it->v0_[0] << " " << it->v0_[1] << " " << it->v0_[2] << std::endl;
    vert << it->v1_[0] << " " << it->v1_[1] << " " << it->v1_[2] << std::endl;
    vert << it->v2_[0] << " " << it->v2_[1] << " " << it->v2_[2] << std::endl;
    face << vnum << ' ' << vnum+1 << ' ' << vnum+2 << std::endl;
    vnum+=3;
  }
}

template <typename T>
int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

Mat3<double> RotationMatrix(double alpha, double beta, double gamma)
{
  double ca = cos(alpha), cb = cos(beta), cg = cos(gamma);
  double sa = sin(alpha), sb = sin(beta), sg = sin(gamma);

  Mat3<double> M;

  M(0,0) = cb*cg;
  M(0,1) = -cb*sg;
  M(0,2) = sb;

  M(1,0) = ca*sg + cg*sa*sb;
  M(1,1) = ca*cg - sa*sb*sg;
  M(1,2) = -cb*sa;

  M(2,0) = sa*sg - ca*cg*sb;
  M(2,1) = cg*sa + ca*sb*sg;
  M(2,2) = ca*cb;

  return M;
}

template <typename VertexType>
void shift(VertexType& v, double xs, double ys, double zs)
{
  v[0] += xs;
  v[1] += ys;
  v[2] += zs;
}

// rotate a point using a rotation matrix
template <typename VertexType, typename MatrixType>
void rotate(VertexType& v, MatrixType& m)
{
  VertexType r(0.);

  r = m*v;
  v = r;
}

template <typename VertexType>
void ConvertRedBloodCellTriangle(VertexType& v)
{
  // constants:
  double r = 3.91, C0 = 0.81, C2 = 7.83, C4 = -4.39;

  double x = v[0]*r, y = v[1]*r;
  double rho = std::sqrt(x*x + y*y);
  double ratio = rho / r;

  // new z co-ordinate
  double z = std::sqrt(1-ratio*ratio + 1e-12)*(C0 + C2*ratio*ratio + C4*ratio*ratio*ratio*ratio)*0.5*sgn(v[2]);

  if (isnan(z)) {
    printf("[E]: z is NaN\n");
    printf("x: %g, y: %g, rho: %g, ratio: %g\n",x,y,rho,ratio);
    double temp = 1-ratio*ratio;
    printf("temp: %g, z: %g\n",temp,z);
    std::exit(0);
  }
  v[0] = x; v[1] = y; v[2] = z;

}

/**
 * Convert a unit sphere into an ethrocyte (RBC)
 */
template <typename PanelType>
void RedBloodCell(std::vector<PanelType>& panels, unsigned recursions = 2, Mat3<double> rotation = Mat3<double>(0.), double *s = nullptr)
{
  std::vector<triangle> triangles(8);
  init_octahedron_triangles(triangles);

  for (unsigned i=0; i<recursions-1; i++) {
    divide_all(triangles);
  }

  // convert all vertices to make RBC
  for (auto it=triangles.begin(); it!=triangles.end(); ++it) {
    // first vertex
    ConvertRedBloodCellTriangle(it->v0_);
    rotate(it->v0_, rotation);
    ConvertRedBloodCellTriangle(it->v1_);
    rotate(it->v1_, rotation);
    ConvertRedBloodCellTriangle(it->v2_);
    rotate(it->v2_, rotation);
    shift(it->v0_, s[0], s[1], s[2]);
    shift(it->v1_, s[0], s[1], s[2]);
    shift(it->v2_, s[0], s[1], s[2]);
  }

  printf("RBC: initialised %d triangles\n",(int)triangles.size());
  // now triangles contains the full triangulation -- convert to panel
  panels.resize(triangles.size());

  auto pit = panels.begin();
  for (auto it=triangles.begin(); it!=triangles.end(); ++it, ++pit) {
    *pit = PanelType(it->v0_,it->v1_,it->v2_);
  }

  // save the triangulation to test.vert, test.face
  std::ofstream face("rbc.face");
  std::ofstream vert("rbc.vert");

  int vnum = 1;
  for (auto it=triangles.begin(); it!=triangles.end();++it) {
    vert << it->v0_[0] << " " << it->v0_[1] << " " << it->v0_[2] << std::endl;
    vert << it->v1_[0] << " " << it->v1_[1] << " " << it->v1_[2] << std::endl;
    vert << it->v2_[0] << " " << it->v2_[1] << " " << it->v2_[2] << std::endl;
    face << vnum << ' ' << vnum+1 << ' ' << vnum+2 << std::endl;
    vnum+=3;
  }
}

/**
 * Create multiple red blood cells, with random orientations
 */
template <typename PanelType>
void MultipleRedBloodCell(std::vector<PanelType>& panels, unsigned recursions = 2, unsigned cells = 1)
{
  unsigned panels_per_cell = 8 * (unsigned)pow(4,recursions-1);
  unsigned total_panels = panels_per_cell * cells;
  panels.resize(total_panels);

  // temp vector to hold a single cell
  std::vector<PanelType> temp;

  unsigned offset = 0;
  double sx = 0., sy = 0., sz = 0.;
  for (unsigned i=0; i<cells; i++) {
    // generate a rotation matrix
    auto M = RotationMatrix(::drand48(), ::drand48(), ::drand48());

    if (i>0) {
      sy += 2*3.91 + 3.91*::drand48();
      sx = i*::drand48() + 2.83*::drand48();
      sz = ((i%2==0)? 1 : -1)*10*::drand48();
    }

    // double sy = 3.91*::drand48() + 2*i*3.91, sx = i*::drand48() + ::drand48(), sz = ((i%2==0)? 1 : -1)*4*::drand48();
    double s[] = {sx, sy, sz};
    // generate the RBC into temp vector
    RedBloodCell(temp, recursions, M, s);

    // shift all panels
    printf("shifting: %g, %g, %g\n",sx,sy,sz);
    /*
    for (auto it=temp.begin(); it!=temp.end(); ++it) {
      shift(it->vertices[0], sx, sy, sz);
      shift(it->vertices[1], sx, sy, sz);
      shift(it->vertices[2], sx, sy, sz);
    }
    */
    // copy temp vector into main panels vector
    for (unsigned j=0; j<panels_per_cell; j++) {
      panels[offset + j] = temp[j];
    }
    // clear the temp vector and increase offset
    temp.clear();
    offset += panels_per_cell;
  }

  // output the entire domain (all cells)
  std::ofstream face("cells.face");
  std::ofstream vert("cells.vert");

  int vnum = 1;
  for (auto it=panels.begin(); it!=panels.end(); ++it)
  {
    auto& v = it->vertices;

    // output vertices
    vert << v[0][0] << " " << v[0][1] << " " << v[0][2] << std::endl;
    vert << v[1][0] << " " << v[1][1] << " " << v[1][2] << std::endl;
    vert << v[2][0] << " " << v[2][1] << " " << v[2][2] << std::endl;
    face << vnum << " " << vnum+1 << " " << vnum+2 << std::endl;
    vnum += 3;
  }
}

}; // end namespace Triangulation
