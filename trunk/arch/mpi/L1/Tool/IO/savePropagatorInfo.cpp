#include <L1/Tool/IO.h>

void Tool::IO::savePropagatorInfo(LemonWriter *writer, size_t L, size_t T)
{
  std::ostringstream oss;
  oss << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n<field>diracFermion</field>\n";
  oss << "<precision>64</precision>\n";
  oss << "<flavours>1</flavours>\n";
  oss << "<lx>" << L << "</lx>\n";
  oss << "<ly>" << L << "</ly>\n";
  oss << "<lz>" << L << "</lz>\n";
  oss << "<lt>" << L << "</lt>\n";
  oss << "</etmcFormat>";
  size_t slen = oss.str().length();

  char *data = new char[slen + 1];
  std::string recordType("etmc-propagator-format");
  std::copy(recordType.begin(), recordType.end(), data);
  LemonRecordHeader *header = lemonCreateHeader(1, 1, data, slen);
  lemonWriteRecordHeader(header, writer);
  lemonDestroyHeader(header);
  std::copy(oss.str().begin(), oss.str().end(), data);
  lemonWriteRecordData(data, &slen, writer);
  delete[] data;
  lemonWriterCloseRecord(writer);
}
