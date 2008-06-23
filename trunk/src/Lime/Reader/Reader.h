#ifndef GUARD_LIME_READER_H
#define GUARD_LIME_READER_H

#include <string>

namespace Lime
{
  struct Data;
  
  class Reader
  {
    Data *d_data;
    int   d_fail;
    
    public:
      Reader(std::string const &filename);
      ~Reader();
      template< typename DataType >
      void read(DataType *buffer, size_t elements) const;
    
    private:
      void read(char *buffer, size_t elements) const;
  };
}

#endif
