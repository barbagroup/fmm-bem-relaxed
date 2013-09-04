#pragma once

/**
 * General purpose mesh reading utilities
 */

namespace MeshIO
{

// split a string by white space
void split(const std::string& str, std::vector<std::string>& v)
{
  std::stringstream ss(str);
  ss >> std::noskipws;

  std::string field;
  char ws_delim;

  while (1) {
    if (ss >> field)
      v.push_back(field);
    else if (ss.eof())
      break;
    else
      v.push_back(std::string());
    ss.clear();
    ss >> ws_delim;
  }
}

}; // end namespace MeshIO
