#include "Reader.ih"

void IO::Lime::Reader::retrieveRecord(int32_t const record)
{
  d_currentRecord = record;
  d_fail = limeSetReaderPointer(d_reader, d_recordOffsets[d_currentRecord]);
}
