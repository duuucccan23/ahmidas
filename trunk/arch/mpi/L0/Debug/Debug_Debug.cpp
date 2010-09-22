/* This output function might seem more complicated then necessary.
 * There is a very good reason for that.
 * Tracing debug information in parallel can be very hard, especially
 * when multiple nodes are writing at the same time. By first building
 * a single string object and dumping that into output in a single
 * call, debug information for individual nodes is preserved.
 * Another way to create such a string would be using the c-style
 * atoi function, but this would reduce flexibility.
 */

#include "Debug.ih"

Debug::Debug(std::string const &debugstring, std::ostream &stream)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ostringstream ostringstream("[DEBUG] Node ", std::ios::ate);
  ostringstream << rank << ": " << debugstring;
  stream << ostringstream.str() << std::endl;
}
