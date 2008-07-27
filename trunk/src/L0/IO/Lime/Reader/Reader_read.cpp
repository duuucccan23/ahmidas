#include "Reader.ih"

void IO::Lime::Reader::read(char *buffer, uint64_t elements) const
{
  if (!elements)
    return; // Lime does not handle this gracefully.
  size_t request = elements;
  d_fail = limeReaderReadData(buffer, &elements, d_reader);
  d_fail = d_fail || (elements != request);
}
