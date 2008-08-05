#include "Writer.ih"

void IO::Lime::Writer::reserveHeader()
{
  for ( size_t ctr = 0; ctr < s_headerSize / 8; ++ctr)
    d_stream.write(s_padding, 8);
  d_record.offset = d_stream.tellp();
}
