#include "Debug.ih"

Debug::Debug(std::string dbgstr)
{
  std::cout << "[DEBUG] " << dbgstr << std::endl << std::flush;
}
