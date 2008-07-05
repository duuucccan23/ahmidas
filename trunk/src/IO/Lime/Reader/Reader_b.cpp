#include "Reader.ih"

void Lime::Reader::read(char *buffer, uint64_t elements) const
{
  if (!elements)
    return; // Lime does not handle this well...
  size_t request = elements;
  d_fail = limeReaderReadData(buffer, &elements, d_data->reader);
  d_fail = d_fail || (elements != request);
}
