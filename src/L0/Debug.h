#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

class Debug
{
  public:
    Debug(std::string const &debugstring, std::ostream &stream = std::cerr);
};
