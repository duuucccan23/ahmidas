#include "Reader.ih"
    
void IO::Lime::Reader::nextRecord()
{
  d_fail = limeReaderNextRecord(d_reader);
}
