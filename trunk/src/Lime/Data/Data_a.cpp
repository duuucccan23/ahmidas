#include "Data.ih"

Lime::Data::Data(std::string const &filename)
  : reader(limeCreateReader(stream))
{
  struct stat statusBuffer;
  int readStatus = stat(filename.c_str(), &statusBuffer);
  
  if (readStatus)
  {
    std::cerr << "An error occurred while attempting to query " << filename << '.' << std::endl;
    MPI::COMM_WORLD.Abort(readStatus);
  }
  
  stream = fopen(filename.c_str(), "r");
  
  if ((!stream) || (!reader))
    MPI::COMM_WORLD.Abort(EIO);
}
