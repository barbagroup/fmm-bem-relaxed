#pragma once

#include <string>
#include <fstream>

namespace DataReaders {

template <typename point_type>
void readVertices(const char *fname, std::vector<point_type>& vertices)
{
  std::ifstream file(fname);
  std::string str;

  while(getline(file, str)) {
    // split the string into x y z
    std::istringstream s(str);
    double x, y, z;
    s >> x >> y >> z;

    // create and append point_type(x,y,z) to array

  }

  file.close();
}

void readTriangle(const char *fname)
{
  std::ifstream file(fname);
  std::string str;

  while (getline(file, str)) {
    // split the string into vertices
    std::istringstream s(str);

    unsigned v1, v2, v3;
    s >> v1 >> v3 >> v2;

    // create and append data structure
  }

  file.close();
}

}; // end namespace DataReaders
