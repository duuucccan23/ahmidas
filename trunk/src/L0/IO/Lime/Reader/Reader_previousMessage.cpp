#include "Reader.ih"
 
void IO::Lime::Reader::previousMessage()
{
  retrieveRecord(d_messageIndices[d_currentRecord] - 1, 0);
}
