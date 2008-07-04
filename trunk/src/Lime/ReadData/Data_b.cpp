#include "Data.ih"

Lime::Data::~Data()
{
  limeDestroyReader(reader);
  fclose(stream);
}
