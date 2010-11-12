#pragma once

#include <stdint.h>
#include <fstream>
#include <string>
#include <cassert>
#include <mpi.h>
#include <L0/Base/Weave.h>

namespace Tool
{
  namespace IO
  {
    namespace Lime
    {
      class Writer
      {

        static size_t const s_headerSize = 144;
        static uint32_t const s_limeMagic = 0x456789ab;
        static char const s_mesBeginMask = 0x80;
        static char const s_mesEndMask = 0x40;
        static char const s_padding[8];

        union uHeader
        {
          uint64_t as64[s_headerSize / sizeof(uint64_t)];
          uint32_t as32[s_headerSize / sizeof(uint32_t)];
          uint16_t as16[s_headerSize / sizeof(uint16_t)];
          char     as8[s_headerSize];
        };

        struct Record
        {
          char           type[128];
          uint64_t       size;
          MPI::Offset    recOffset; // absolute offset at which the record starts
          MPI::Offset    offset;
          int16_t        version;
          bool           mesBeg;
          bool           mesEnd;

          Record();
          Record(size_t const rOffset, size_t const recSize = 0);
        };

        private:

          Base::Weave    *d_weave;
          MPI::File      d_MPI_FILE;
          MPI::Offset    d_startOfNextRecord;
          MPI::Status    d_MPI_Status;
          Record         d_record;
          bool           d_hasWritten;
          bool           d_messageRunning;
          bool           d_writeHeader;

        public:
          Writer(Base::Weave passed_weave,std::string const &filename);
          ~Writer();

          void finishMessage();

          // for parallel writing also the size should be passed. If not (or zero is passed)
          // the function assumes that only onl process writes the record
          void newRecord(std::string const &type, size_t const rOffset, size_t const size = 0);

          size_t closeRecord();

          template< typename DataType >
          void write(DataType const *buffer, uint64_t elements);

          template< typename DataType >
          void write(DataType const *buffer, DataType const *finish);

          void write_collective(double const *buffer, uint64_t const count, uint64_t const double_offset);

          void preallocate(uint64_t bytes);
          void reset_view();

          bool fail() const;
          bool good() const;

          // sets position to offset bytes in current record

          //void seekp(MPI::Offset const offset);
          MPI::Offset tellp();

        private:
          void finalize();
          void reserveHeader();
      };
    }
  }
}

#include "Writer/Writer.inlines"
