#include "WriteData.ih"

Lime::WriteData::~WriteData()
{
  limeDestroyWriter(writer);
  fclose(stream);
}
