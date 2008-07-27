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
    if (d_reader.limeType() == std::string("ildg-data-lfn")) // Compare differently?
    {
      d_lfnRecord = d_reader->currentRecord();
      char *lfn = new char[d_reader.size()];
      d_reader->read(lfn, d_reader->size());
      d_lfn = std::string(lfn);
      delete lfn;
      continue;
    }
    if (d_reader.limeType() == std::string("ildg-format"))
  }

  return d_reader.fail();
}
