#pragma once

#include <string>

namespace IO
{
  struct ILDGFormat
  {
    std::string ildgFormat;
    std::string version;
    std::string field;
  
    size_t nx;
    size_t ny;
    size_t nz;
    size_t nt;
    size_t precision;
    
    void parse(char *message);
  };
}
