#include "ReadData.ih"

Lime::ReadData::~ReadData()
{
  limeDestroyReader(reader);
  fclose(stream);
}
