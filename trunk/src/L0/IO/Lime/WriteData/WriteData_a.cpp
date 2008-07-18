#include "WriteData.ih"

Lime::WriteData::WriteData(std::string const &filename)
{
  stream = fopen(filename.c_str(), "w");
  
 if (!stream)
    MPI::COMM_WORLD.Abort(EIO);
  
  writer = limeCreateWriter(stream);

  if (!writer)
    MPI::COMM_WORLD.Abort(EIO);
}
