#include "Print.ih"

Print::Print(std::ostream *strm, std::string printstr)
{
  if (!MPI::COMM_WORLD.Get_rank())
    *strm << printstr << std::endl;
}
