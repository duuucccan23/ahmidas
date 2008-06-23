#ifndef GUARD_LIME_READER_H
#define GUARD_LIME_READER_H

#include <string>

namespace Lime
{
  struct Data;
  
  class Reader
  {
    Data        *d_data;
    int mutable  d_fail;
    size_t       d_size;
    
    public:
      Reader(std::string const &filename);
      ~Reader();
      
      int fail() const;
      size_t size() const;
      
      template< typename DataType >
      void read(DataType *buffer, size_t elements) const;
    
    private:
      void read(char *buffer, size_t elements) const;
  };
}

#include "Reader.inlines"

#endif
