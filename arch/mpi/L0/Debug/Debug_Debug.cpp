#include "Debug.ih"

Debug::Debug(std::string dbgstr)
{
  std::cout << "[DEBUG] Node " << MPI::COMM_WORLD.Get_rank() << ": " << dbgstr << std::endl;
}
