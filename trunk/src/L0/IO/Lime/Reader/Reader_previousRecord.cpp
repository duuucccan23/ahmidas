#include "Reader.ih"

void IO::Lime::Reader::previousRecord()
{
  d_fail = limeSetReaderPointer(d_reader, d_recordOffsets[--d_currentRecord]);
}
