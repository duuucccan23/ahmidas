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
  oss << "<lt>" << T << "</lt>\n";
  oss << "</etmcFormat>";
  std::string const datastr(oss.str());
  size_t slen = datastr.length();

  char *data = new char[slen + 1];
  std::fill_n(data, slen+1, '\0');
  std::string recordType("etmc-propagator-format");
  std::copy(recordType.begin(), recordType.end(), data);
  LemonRecordHeader *header = lemonCreateHeader(1, 1, data, slen);
  lemonWriteRecordHeader(header, writer);
  lemonDestroyHeader(header);
  std::copy(datastr.begin(), datastr.end(), data);
  lemonWriteRecordData(data, &slen, writer);
  delete[] data;
  lemonWriterCloseRecord(writer);
}
