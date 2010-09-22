#include "Print.ih"

Print::Print(std::string const &printstr, std::ostream &strm)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    strm << printstr << std::endl;
}
