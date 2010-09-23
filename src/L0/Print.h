#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

class Print
{
  public:
    Print(std::string const &printstr, std::ostream &strm = std::cout);
};
