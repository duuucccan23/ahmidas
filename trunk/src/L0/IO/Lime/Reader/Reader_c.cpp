#include "Reader.ih"

Lime::Reader::~Reader()
{
  limeDestroyWriter(d_reader);
  fclose(d_stream);
}
