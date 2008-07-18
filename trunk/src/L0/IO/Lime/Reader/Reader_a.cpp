#include "Reader.ih"

Lime::Reader::Reader(std::string const &filename, std::string const &filetype)
  : d_data(new ReadData(filename)), d_fail(0), d_size(0)
{
  while(!(d_fail = limeReaderNextRecord(d_data->reader)))
  {
    d_data->headerType = limeReaderType(d_data->reader);
    if(!strcmp(filetype.c_str(), d_data->headerType))
      break;
  }

  switch(d_fail)
  {
    case LIME_SUCCESS:
      break;
    case LIME_EOF:
      std::cerr << "The requested file " << filename << " was not recognized as being of type " << filetype << std::endl;
      MPI::COMM_WORLD.Abort(EIO);
    default:
      std::cerr << "Lime reader exited with I/O error " << d_fail << '.' << std::endl;
      MPI::COMM_WORLD.Abort(EIO);
  }

  d_size = limeReaderBytes(d_data->reader);
}
