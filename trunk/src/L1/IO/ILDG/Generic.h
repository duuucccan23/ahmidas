#ifndef GUARD_IO_ILDG_GENERIC
#define GUARD_IO_ILDG_GENERIC

#include <string>

#include <L0/Core/Core.h>
#include <L0/IO/Lime/Reader.h>

namespace IO
{
  namespace ILDG
  {
    struct Format
    {
      std::string version;
      std::string field;
      size_t      precision;
      size_t      dimensions[4];
    };

    class Generic
    {
      std::string const  d_filename;
      Lime::Reader      *d_reader;

      Format             d_format;
      std::string        d_lfn;

      uint32_t           d_payloadMessage;
      uint32_t           d_formatRecord;
      uint32_t           d_dataRecord;
      uint32_t           d_lfnRecord;

      public:
        Generic(std::string const &filename);
        void openFile(std::string const &filename);
        int fail() const;

        template< typename DataType >
        void read(DataType *buffer, uint64_t const elements) const;
        
        std::string const &version() const;
        std::string const &field() const;
        size_t const precision() const;
        size_t const *dimensions() const;

      private:
        void parseXmlFormat(char *data);
    };
  }
}

#include "Generic/Generic.inlines"
#include "Generic/Generic_read.template"

// Add a specialization for the Core::Field::readDataFromIO method if needed
#ifdef FLAG_L0_CORE_FIELD_IS_KNOWN 
#include "Generic/Field_loadDataFromIO.template"
#endif
             
#define FLAG_L1_IO_ILDG_GENERIC_IS_KNOWN
 
#endif
