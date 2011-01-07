#include <L1/Tool/IO.h>

void Tool::IO::savePropagatorType(LemonWriter *writer)
{
  std::string infoString("DiracFermion_Sink");
  uint64_t slen = infoString.length();

  char data[32];
  std::string recordType("propagator-type");
  std::copy(recordType.begin(), recordType.end(), data);
  LemonRecordHeader *header = lemonCreateHeader(1, 1, data, slen);
  lemonWriteRecordHeader(header, writer);
  lemonDestroyHeader(header);
  std::copy(infoString.begin(), infoString.end(), data);
  lemonWriteRecordData(data, &slen, writer);
  lemonWriterCloseRecord(writer);
}
