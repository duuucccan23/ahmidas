#include "Debug.ih"

Debug::Debug(std::string const &debugstring, std::ostream &stream)
{
  stream << "[DEBUG] : " << debugstring << std::endl;
}
