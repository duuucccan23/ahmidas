#include <L1/Tool/IO.h>

void Tool::IO::loadScidacInfo(LemonReader *reader, ScidacInfo &info)
{
  uint64_t nbytes = lemonReaderBytes(reader);
  char *scidacCStr = new char[nbytes];
  lemonReaderReadData(scidacCStr, &nbytes, reader);
  info = parseScidacInfo(scidacCStr);
  delete[] scidacCStr;
}
