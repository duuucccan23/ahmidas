#include "Print.ih"

Print::Print(std::string const &printstr, std::ostream &strm)
{
  strm << printstr << std::endl;
}
