#include "Print.ih"

Print::Print(std::ostream *strm, std::string printstr)
{
  *strm << printstr << std::endl;
}
