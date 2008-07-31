#include "Generic.ih"

void IO::ILDG::Generic::openFile(std::string const &filename)
{
  delete d_reader;
  d_reader = new Lime::Reader(filename);

  std::vector< std::string > const &types = d_reader->limeTypes();
  d_formatRecord = std::find(types.begin(), types.end(), std::string("ildg-format")) - types.begin();
  d_dataRecord = std::find(types.begin(), types.end(), std::string("ildg-binary-data")) - types.begin();
  d_lfnRecord = std::find(types.begin(), types.end(), std::string("ildg-data-lfn")) - types.begin();
  
  if (d_formatRecord == types.size() || d_dataRecord == types.size())
    MPI::COMM_WORLD.Abort(EIO); // NOTE Invalid ILDG file, mandatory fields missing
  
  d_reader->retrieveRecord(d_formatRecord);
  
  d_payloadMessage = d_reader->currentMessage();
  std::cerr  << "[DEBUG] Record:" << d_reader->currentRecord() << std::endl;
    
  char *format = new char[d_reader->size() + 1];
  d_reader->read(format, d_reader->size());
  format[d_reader->size()] = '\0';
  std::cerr << "[DEBUG] setup:" << d_reader->size() << "  " << format << std::endl;
  parseXmlFormat(format);
  delete[] format;

  if (d_lfnRecord == types.size())
  {
    std::cerr << "Mandatory LFN LIME record missing from file!\n"
              << "See ILDG Binary File Format description rev 1.1 for details.\n"
              << "Permitting file read in for now, however." << std::endl;
  }
  else
  {
    d_reader->retrieveRecord(d_lfnRecord);
    char *lfn = new char[d_reader->size()];
    d_reader->read(lfn, d_reader->size());
    d_lfn = std::string(lfn);
    delete[] lfn;
  }
  
  d_reader->retrieveRecord(d_dataRecord); // Be kind, rewind :)
}
