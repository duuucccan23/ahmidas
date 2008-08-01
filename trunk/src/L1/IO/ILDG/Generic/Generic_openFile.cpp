#include "Generic.ih"

void IO::ILDG::Generic::openFile(std::string const &filename)
{
  delete d_reader;
  d_reader = new Lime::Reader(filename);

  d_formatRecord = d_reader->findRecord("ildg-format");
  d_dataRecord = d_reader->findRecord("ildg-binary-data");
  d_lfnRecord = d_reader->findRecord("ildg-data-lfn");
  
  if (d_formatRecord == d_reader->records() || d_dataRecord == d_reader->records())
    MPI::COMM_WORLD.Abort(EIO); // NOTE Invalid ILDG file, mandatory fields missing
  
  d_reader->retrieveRecord(d_formatRecord);
  
  d_payloadMessage = d_reader->currentMessage();
  std::cerr  << "[DEBUG] Record:" << d_reader->currentRecord() << std::endl;
    
  char *format = new char[d_reader->recordSize() + 1];
  d_reader->read(format, d_reader->recordSize());
  format[d_reader->recordSize()] = '\0';
  std::cerr << "[DEBUG] setup:" << d_reader->recordSize() << "  " << format << std::endl;
  parseXmlFormat(format);
  delete[] format;

  if (d_lfnRecord == d_reader->records())
  {
    std::cerr << "Mandatory LFN LIME record missing from file!\n"
              << "See ILDG Binary File Format description rev 1.1 for details.\n"
              << "Permitting file read in for now, however." << std::endl;
  }
  else
  {
    d_reader->retrieveRecord(d_lfnRecord);
    char *lfn = new char[d_reader->recordSize()];
    d_reader->read(lfn, d_reader->recordSize());
    d_lfn = std::string(lfn);
    delete[] lfn;
  }
  
  d_reader->retrieveRecord(d_dataRecord); // Be kind, rewind :)
}
