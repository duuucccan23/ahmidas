#include "Data.ih"

Lime::Data::Data(std::string const &filename)
  : reader(limeCreateReader(stream))
{
  std::cerr << "[DEBUG] Entered Lime::Data constructor" << std::endl;
  struct stat statusBuffer;
  int readStatus = stat(filename.c_str(), &statusBuffer);
  std::cerr << "[DEBUG] Performed stat on " << filename << std::endl;
  
  switch (readStatus)
  {
    case 0:
      std::cerr << "[DEBUG] Managed to read file normally" << std::endl;
      break;
    case EACCES:
      std::cerr << "An error occured while trying to access " << filename 
		<< ", make sure it has read permission!" << std::endl;
      break;
    case ENOENT:
      std::cerr << "It seems " << filename << " does not exist." << std::endl;
      break;
    case ENAMETOOLONG:
      std::cerr << "The name " << filename << " is too long to be a proper filename." << std::endl;
      break;
    case ENOMEM:
      std::cerr << "Internal error: kernel memory depleted when attempting to query " << filename
		<< '.' << std::endl;
      break;
    default:
      std::cerr << "An unexpected error occurred while attempting to query " 
		<< filename << '.' << std::endl;
      break;
  }
  
  std::cerr << "[DEBUG] Preparing to set up stream" << std::endl;
  
  if (readStatus)
    MPI::COMM_WORLD.Abort(readStatus);
  
  stream = fopen(filename.c_str(), "r");
  
  if ((!stream) || (!reader))
    MPI::COMM_WORLD.Abort(EIO);
}
