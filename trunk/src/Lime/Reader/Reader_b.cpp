#include "Reader.ih"

void Lime::Reader::read(double *buffer, size_t elements) const
{
  size_t request = elements;
  d_fail = limeReaderReadData(buffer, &elements, d_data->reader);
  d_fail = d_fail && (elements == request);
}
