#pragma once

#include <string>
#include <iostream>

class Debug
{
  public:
    Debug(std::string const &debugstring, std::ostream &stream = std::cerr);
};
