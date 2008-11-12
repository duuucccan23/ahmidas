#include "MILCinfo.ih"

Tool::IO::MILCinfo::MILCinfo(std::ifstream &reader)
{
  static_assert(sizeof(Tool::IO::MILCinfo) == 96, "MILCinfo stored in an unexpected way.");
  reader.read(reinterpret_cast< char* >(this), 96);
  if (!Base::bigEndian)
    Base::swapEndian(&checksum, &checksum + 1);
}
