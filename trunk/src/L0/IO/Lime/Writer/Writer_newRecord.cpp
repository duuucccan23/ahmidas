#include "Writer.ih"

void IO::Lime::Writer::newRecord(std::string const &type)
{
  finalize();

  bool startMessage = d_record.mesEnd;
  d_record = Record();
  d_record.mesBeg = startMessage;

  size_t typeSize = type.size() < 127 ? type.size() : 127;
  std::copy(type.c_str(), type.c_str() + typeSize, d_record.type);

  reserveHeader();
}
