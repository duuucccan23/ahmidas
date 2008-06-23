#include "Data.ih"

Lime::Data::~Data()
{
  limeDestroyReader(d_data->reader);
  fclose(d_data->stream);
}
