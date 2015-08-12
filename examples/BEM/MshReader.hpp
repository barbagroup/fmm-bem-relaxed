#pragma once

/**
 * Read in a gmsh .msh format and output as vector of triangles
 */

#include <vector>
#include <cstring>
#include <string>
#include <fstream>

#include "MeshIO.hpp"

namespace MeshIO
{


template <typename point_type, typename triangle_type>
void readMsh(const char *fname, std::vector<triangle_type>& elements)
{
  std::ifstream mesh(fname);

  // read header line
  std::string str;
  // information about the mesh
  int num_nodes = 0, num_elements = 0;
  std::vector<point_type> nodes;

  bool in_header = true;
  while (in_header) {
    std::getline(mesh, str);
    if (str == "$Nodes") in_header = false;
  }

  // this is # of nodes
  std::getline(mesh,str);
  num_nodes = atoi(str.c_str());
  nodes.resize(num_nodes);
  printf("num_nodes: %d\n",num_nodes);

  std::vector<std::string> split_str;
  // now read all num_nodes nodes into an array
  for (int i=0; i<num_nodes; i++) {
    split_str.clear();
    std::getline(mesh, str); // read current line

    split(str.c_str(), split_str);

    // grab data from this node
    int node_no;
    double x, y, z;
    node_no = atoi(split_str[0].c_str());
    x = atof(split_str[1].c_str());
    y = atof(split_str[2].c_str());
    z = atof(split_str[3].c_str());
    // create a node
    nodes[node_no-1] = point_type(x,y,z); // adjust for 1-indexing
  }
  // now read $EndNodes and $Elements
  std::getline(mesh, str);
  std::getline(mesh, str);
  // this line is # elements
  std::getline(mesh, str);
  num_elements = atoi(str.c_str());
  printf("num_elements: %d\n",num_elements);
  elements.resize(num_elements);

  int skipped = 0; // # skipped elements
  for (int i=0; i<num_elements; i++) {
    split_str.clear();

    std::getline(mesh, str);
    split(str.c_str(), split_str);

    int element_no = atoi(split_str[0].c_str());
    if (atoi(split_str[1].c_str()) != 2) {
      skipped++;
      // printf("\tNon-triangular element found -- skipping\n");
      continue;
    }
    // int element_type = atoi(split_str[1].c_str());
    int num_properties = atoi(split_str[2].c_str());

    int v1 = atoi(split_str[2+num_properties+1].c_str()) - 1; // adjust for 1-indexing
    int v2 = atoi(split_str[2+num_properties+2].c_str()) - 1;
    int v3 = atoi(split_str[2+num_properties+3].c_str()) - 1;

    // elements[element_no-1] = triangle_type(nodes[v1],nodes[v2],nodes[v3]); // clockwise
    elements[element_no-1] = triangle_type(nodes[v1],nodes[v3],nodes[v2]); // anti-clockwise
  }

  printf("%d elements skipped\n",skipped);
  elements.resize(num_elements-skipped);
}

}; // end namespace MeshIO
