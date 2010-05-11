#include "Writer.ih"
size_t Tool::IO::Lime::Writer::closeRecord()
{
  if (!d_hasWritten)
    finalize();
  d_hasWritten = true;

  if (d_writeHeader)
  {
    // std::cout << "d_startOfNextRecord = " << d_startOfNextRecord << std::endl;
    return d_startOfNextRecord;
  }
  else return 0;
}
