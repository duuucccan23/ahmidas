#include "Reader.ih"

void Lime::Reader::read(double *buffer, size_t elements) const
{
  d_fail = limeReaderReadData(buffer, elements, d_data->reader);
}
