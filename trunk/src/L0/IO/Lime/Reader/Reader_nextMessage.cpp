#include "Reader.ih"

void IO::Lime::Reader::nextMessage()
{
  retrieveRecord(d_messageIndices[d_currentRecord] + 1, 0);
}
