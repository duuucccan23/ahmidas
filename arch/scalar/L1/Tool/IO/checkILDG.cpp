#include "../IO.h"

int Tool::IO::checkILDG(std::string const &filename)
{
  Lime::Reader reader(filename);
  if(!reader.good())
  {
    std::cerr << "Error occured trying to open file: " << filename << std::endl;
    return (-2);
  }
  assert(reader.good());

  ILDGinfo info(reader);
  reader.retrieveRecord(reader.findRecord("ildg-binary-data"));
  if (reader.fail())
  {
    std::cerr << "Lime reader could not find a ildg-binary-data record in file: " << filename << std::endl;
    return (-3);
  }
  assert(reader.good());

  QCD::Gauge buf;
  Tool::ScidacChecksum calc(0);

  size_t siteItr = 0;
  size_t sites = info.dims[0] * info.dims[1] * info.dims[2] * info.dims[3];
  while (siteItr < sites) {
    reader.read(&buf, 1);
    calc.aggregate(&buf, 1, siteItr);
    ++siteItr;
  }

  reader.retrieveRecord(reader.findRecord("scidac-checksum"));
  if (reader.fail())
  {
    std::cerr << "Lime reader could not find a scidac-checksum record in file: " << filename << std::endl;
    return (-4);
  }
  assert(reader.good());

  char mymesg[reader.recordSize()];
  reader.read(mymesg, reader.recordSize());

  uint32_t reada;
  uint32_t readb;
  char *pos = strtok(mymesg, "<> \n\t");

  while (pos)
  {
    if (!strncmp(pos, "suma", 4)) {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%x", &reada);
    }
    if (!strncmp(pos, "sumb", 4)) {
      pos = strtok(0, "<> \n\t");
      sscanf(pos, "%x", &readb);
    }
    pos = strtok(0, "<> \n\t");
  }
  if ((reada != calc.lower()) | (readb != calc.upper()))
    return(-1);
  return(0);
}
