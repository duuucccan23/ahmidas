#include "Generic.ih"

void IO::ILDG::Generic::readFile(std::string const &filename)
{
  delete d_reader;
  d_reader = new Lime::Reader(filename);
  if (d_reader->messages() < 2)
    MPI::COMM_WORLD.Abort(EIO);

  while (!d_reader->fail())
  {
    d_reader->nextRecord();
    if (d_reader.limeType() == std::string("ildg-binary-data"))
    {
      d_dataRecord = d_reader->currentRecord();
      continue;
    }
    if (d_reader.limeType() == std::string("ildg-data-lfn"))
    {
      d_lfnRecord = d_reader->currentRecord();
      char *lfn = new char[d_reader.size()];
      d_reader->read(lfn, d_reader->size());
      d_lfn = std::string(lfn);
      delete[] lfn;
      continue;
    }
    if (d_reader.limeType() == std::string("ildg-format"))
    {
      d_payloadMessage = d_reader->currentMessage();
      d_formatRecord = d_reader->currentRecord();
      char *format = new char[d_reader.size()];
      d_reader->read(format, d_reader->size()];
      parseXmlFormat(format);
      delete[] format;
    }
  }
  d_reader->retrieveRecord(d_dataRecord); // Be kind, rewind :)
}
