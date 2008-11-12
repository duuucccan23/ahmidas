#include "MILCinfo.ih"

Tool::IO::MILCinfo::MILCinfo(std::ifstream &reader)
{
  static_assert(sizeof(Tool::IO::MILCinfo) == 96, "MILCinfo stored in an unexpected way.");
  reader.read(reinterpret_cast< char* >(this), 96);
  std::cerr << "[DEBUG] Contents of MILCinfo after parsing" << std::endl;
  std::cerr << "\tmagicNumber\t" << magicNumber << std::endl;
  std::cerr << "\ttimestamp\t" << timestamp << std::endl;
  std::cerr << "\tdims\t" << dims[0] << ' ' << dims[1] << ' ' << dims[2] << ' ' << dims[3] << std::endl;
  std::cerr << "\torder\t" << order << std::endl;
  std::cerr << "\tchecksum\t" << checksum << std::endl;
}
