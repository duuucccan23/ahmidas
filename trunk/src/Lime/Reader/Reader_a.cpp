#include "Reader.ih"

void Lime::Reader(std::string const &filename) const
  : d_data(new Data(filename)), d_fail(0), d_size(0)
{
  while(!(d_fail = limeReaderNextRecord(d_data->reader)))
  {
    d_data->headerType = limeReaderType(limereader);
    if(!strcmp("ildg-binary-data", header_type))
      break;
  }
  
  swith(d_fail)
  {
    case LIME_SUCCESS:
      break;
    case LIME_EOF:
      std::cerr << "No ildg-binary-data record found in " << filename << '.' << endl;
      MPI::COMM_WORLD.Abort(EIO);
    default:
      std::cerr << "Lime reader exited with I/O error " << d_fail << '.' << endl;
      MPI::COMM_WORLD.Abort(EIO);
  }
  
  d_size = limeReaderBytes(d_data->reader);
}
