#include "Writer.ih"

void Tool::IO::Lime::Writer::reserveHeader()
{
  if(!d_writeHeader)
    return;

  // d_record.offset should be zero, and position in file should be d_record.recOffset
  assert(d_stream.tellp() == d_record.recOffset);
  assert(d_record.offset == std::streampos(0));

  for ( size_t ctr = 0; ctr < s_headerSize / 8; ++ctr)
    d_stream.write(s_padding, 8);
}
