#include "Reader.ih"

IO::Lime::Reader::~Reader()
{
  limeDestroyReader(d_reader);
  fclose(d_stream);
}
