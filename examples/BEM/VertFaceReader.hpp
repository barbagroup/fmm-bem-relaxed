#pragma once

/**
 * Read in a mesh defined as .vert and .face
 */

#include <vector>
#include <cstring>
#include <string>
#include <fstream>

#include "MeshIO.hpp"

namespace MeshIO
{

template <typename point_type, typename triangle_type>
void ReadVertFace(const char *vert_file, const char *face_file, std::vector<triangle_type>& elements)
{
  std::ifstream vert(vert_file);
  std::ifstream face(face_file);

  std::string str;

  // mesh info
  int num_vertices = 0, num_faces = 0;
  std::vector<point_type> vertices;

  // get first line of vertices -- # vertices
  std::getline(vert, str);
  num_vertices = atoi(str.c_str());
  vertices.resize(num_vertices);
  printf("# vertices: %d\n",num_vertices);

  std::vector<std::string> split_str;
  // read all vertices
  for (int i=0; i<num_vertices; i++) {
    split_str.clear();
    std::getline(vert, str);
    split(str.c_str(), split_str);

    double x, y, z;
    x = atof(split_str[0].c_str());
    y = atof(split_str[1].c_str());
    z = atof(split_str[2].c_str());

    vertices[i] = point_type(x,y,z);
  }
  // finished reading vertices
  vert.close();

  // get first line of faces -- # faces
  std::getline(face, str);
  num_faces = atoi(str.c_str());
  elements.resize(num_faces);
  printf("# elements: %d\n",num_faces);

  // read all faces
  for (int i=0; i<num_faces; i++) {
    split_str.clear();
    std::getline(face, str);
    split(str.c_str(), split_str);

    int v1, v2, v3;
    // read vertex numbers, taking into account 1-indexing
    v1 = atoi(split_str[0].c_str()) - 1;
    v2 = atoi(split_str[1].c_str()) - 1;
    v3 = atoi(split_str[2].c_str()) - 1;

    elements[i] = triangle_type(vertices[v1],vertices[v2],vertices[v3]);
  }
  face.close();
}

}; // end namespace MeshIO

