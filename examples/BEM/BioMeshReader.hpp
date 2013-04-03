#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

/** Class to read a mesh from following files:
 *  .vert |  x  y  z <not needed>
 *  .face | v1 v3 v2 <not needed>
 * into panels */
template <class PointType, class PanelType>
class BioMeshReader
{
 public:
  static void ReadMesh(const char *vert_file, const char *face_file, std::vector<PanelType>& panels) {

    // temp storage for vertices
    std::vector<PointType> vertices;

    std::ifstream vert(vert_file);
    std::ifstream face(face_file);

    std::string line;
    // read vertices
    while (getline(vert,line)) {
      std::istringstream i(line);

      double x, y, z;
      i >> x >> y >> z;

      // create panel -- for now just print
      vertices.push_back(PointType(x,y,z));
    }
    vert.close();

    while (getline(face,line)) {
      std::istringstream i(line);

      int v1, v2, v3;
      i >> v1 >> v3 >> v2;

      // deal with 0 vs. 1-indexing
      v1--; v2--; v3--;

      // append this panel
      panels.push_back(PanelType(vertices[v1], vertices[v2], vertices[v3]));
    }

    face.close();
  }

};
