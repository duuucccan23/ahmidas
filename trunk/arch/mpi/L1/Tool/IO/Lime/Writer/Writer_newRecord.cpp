#include "Writer.ih"

void Tool::IO::Lime::Writer::newRecord(std::string const &type, size_t const rOffset, size_t const size)
{

  d_record = Record(rOffset, size);
  if (!d_messageRunning)
    d_record.mesBeg = d_messageRunning = true;

  size_t typeSize = type.size() < 127 ? type.size() : 127;
  std::copy(type.c_str(), type.c_str() + typeSize, d_record.type);

  d_startOfNextRecord = rOffset;

  if (d_writeHeader)
   d_MPI_FILE.Seek(d_startOfNextRecord, MPI_SEEK_SET);

  reserveHeader();
}
