#ifndef GUARD_LIME_WRITER_H
#define GUARD_LIME_WRITER_H

namespace Lime
{
  class Writer
  {
    WriteData   *d_data;
    int mutable  d_fail;
    uint64_t     d_written;
    uint64_t     d_totalSize;

    public:
      Writer(std::string const &filename, std::string const &filetype);
      ~Writer();

      int fail() const;
      size_t size() const;

      template< typename DataType >
      void setupFile(uint64_t ) const
      
      void write(char *buffer, uint64_t elements) const;

      template< typename DataType >
      void write(DataType *buffer, uint64_t elements) const;    
  };
}

#include "Writer.inlines"

#endif
