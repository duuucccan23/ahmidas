#include "Reader.ih"

void IO::Lime::Reader::retrieveRecord(int32_t const message, int32_t const relativeRecord)
{
  int32_t idx = std::find(d_messageIndices.begin(), d_messageIndices.end(), message) - d_messageIndices.begin();
  d_currentRecord = idx + relativeRecord;
  d_fail = limeSetReaderPointer(d_reader, d_recordOffsets[d_currentRecord]);
}
