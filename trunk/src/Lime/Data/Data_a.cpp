#include "Data.ih"

Lime::Data::Data(std::string const &filename)
  : stream(fopen(filename.c_str(), "r")), reader(limeCreateReader(stream))
{
  if ((!stream) || (!reader))
    MPI::COMM_WORLD.Abort(EIO);
}
