#include "MILCinfo.ih"

Tool::IO::MILCinfo::MILCinfo(std::ifstream &reader)
{
  reader.read(reinterpret_cast< char* >(this), 96);
  if (!Base::bigEndian)
    Base::swapEndian(&checksum, &checksum + 1);
}
