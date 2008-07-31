#ifndef GUARD_IO_GENERIC_H
#define GUARD_IO_GENERIC_H

#include <fstream>
#include <string>

namespace IO
{
  class Generic
  {
    std::ifstream d_stream;
    
    public:
      Generic(std::string const &filename);
      
      template< typename DataType >
      void read(DataType *buffer, uint64_t const elements) const;
      
      bool fail() const;
  };
}

#include "Generic/Generic.inlines"

#endif
