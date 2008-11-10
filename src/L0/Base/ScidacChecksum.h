#ifndef GUARD_BASE_SCIDACCHECKSUM_H
#define GUARD_BASE_SCIDACCHECKSUM_H

#include <cstdint>
#include <cstdlib>

namespace Base
{
  class ScidacChecksum
  {
    union
    {
      uint32_t as32[2];
      uint64_t as64;
    } d_sum;

    uint32_t d_crcTable[256];

    public:
      ScidacChecksum();
      void clear();
      uint64_t checksum();

      template< typename Element >
      size_t aggregate(Element *data, size_t elements, size_t rank = 0);

    private:
      void createTable();

      template< typename Element >
      uint32_t crc32(Element &buffer);

      template< typename Element >
      uint32_t crc32(Element *buffer, size_t length = 1);
  };
}

#include "ScidacChecksum/ScidacChecksum.inlines"
#include "ScidacChecksum/ScidacChecksum_aggregate.template"
#include "ScidacChecksum/ScidacChecksum_crc32_a.template"
#include "ScidacChecksum/ScidacChecksum_crc32_b.template"

#endif
