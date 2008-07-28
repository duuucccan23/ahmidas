#include "Generic.ih"

void IO::ILDG::Generic::parseXmlFormat(char *data)
{
  std::istringstream xmlFormat(data);
  std::string word;
  while (true)
  {
    xmlFormat >> word;
    if (!xmlFormat)
      return; // String has been fully read.
    if (word == std::string("<version>"))
    {
      xmlFormat >> d_format.version;
      continue;
    }
    if (word == std::string("<field>"))
    {
      xmlFormat >> d_format.field;
      continue;
    }
    if (word == std::string("<precision>"))
    {
      xmlFormat >> d_format.precision;
      continue;
    }
    if (word == std::string("<lx>"))
    {
      xmlFormat >> d_format.dimensions[Core::idx_X];
      continue;
    }
    if (word == std::string("<ly>"))
    {
      xmlFormat >> d_format.dimensions[Core::idx_Y];
      continue;
    }
    if (word == std::string("<lz>"))
    {
      xmlFormat >> d_format.dimensions[Core::idx_Z];
      continue;
    }
    if (word == std::string("<lt>"))
      xmlFormat >> d_format.dimensions[Core::idx_T];
  }
}
