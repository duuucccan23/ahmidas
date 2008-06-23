#include "Reader.ih"

Lime::Reader::Reader(std::string const &filename)
  : d_data(new Data(filename)), d_fail(0), d_size(0)
{
  while(!(d_fail = limeReaderNextRecord(d_data->reader)))
  {
    d_data->headerType = limeReaderType(d_data->reader);
    if(!strcmp("ildg-binary-data", d_data->headerType))
      break;
  }
  
  switch(d_fail)
  {
    case LIME_SUCCESS:
      break;
    case LIME_EOF:
      std::cerr << "No ildg-binary-data record found in " << filename << '.' << std::endl;
      MPI::COMM_WORLD.Abort(EIO);
    default:
      std::cerr << "Lime reader exited with I/O error " << d_fail << '.' << std::endl;
      MPI::COMM_WORLD.Abort(EIO);
  }
  
  d_size = limeReaderBytes(d_data->reader);
}
