#include "Reader.ih"

IO::Lime::Reader::Reader(std::string const &filename)
  : d_filename(filename), d_in(filename.c_str()), d_currentRecord(0),
    d_fail(false), d_mesTerm(true)
{
  uHeader header;  
  bool bigEndian = Core::big_endian();
  int32_t message = -1;
  while (d_in.good())
  {
    d_in.read(header.as8, s_headerSize);
    
    if (!bigEndian)
      Core::swapEndian(header.as32[0]);
    if (header.as32[0] != s_limeMagic)
      break; // We're done, apparently.
        
    if (!bigEndian)
      Core::swapEndian(header.as16[2]);
    d_versions.push_back(header.as16[2]);
    
    if (d_mesTerm && (header.as8[6] & s_mesBeginMask))
    {
      ++message;
      d_mesTerm = false;
    }
    if (!d_mesTerm && (header.as8[6] & s_mesEndMask))
      d_mesTerm = true;
    d_messages.push_back(message);
    
    if (!bigEndian)
      Core::swapEndian(header.as64[1]);
    d_sizes.push_back(header.as64[1]);
    
    d_types.push_back(std::string(&header.as8[16]));
        
    d_offsets.push_back(d_in.tellg());   
    d_in.seekg((d_sizes.back() + 7) / 8 * 8 , std::ios::cur);
  }
  d_in.clear();
  d_in.seekg(d_offsets.front(), std::ios::beg);
}
