#pragma once

#include <string>
#include <iostream>

class Print
{
  public:
    Print(std::string const &printstr, std::ostream &strm = std::cout);
};
