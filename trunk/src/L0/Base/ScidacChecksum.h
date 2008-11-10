#ifndef GUARD_BASE_SCIDACCHECKSUM_H
#define GUARD_BASE_SCIDACCHECKSUM_H

#include <L0/Core/Field.h>

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

      uint64_t checksum() const;
      uint32_t lower() const;
      uint32_t upper() const;

      template< typename Element >
      size_t aggregate(Element const *data, size_t elements, size_t rank = 0);

      template< typename Element >
      size_t blockAggregate(Element const *data, size_t blockSize, size_t rank = 0, size_t blocks = 1);

      template< typename Element, size_t L, size_t T >
      void calculate(Core::Field< Element, L, T > const &field);

    public:
      void createTable();

      template< typename Element >
      uint32_t crc32(Element const &buffer);

      template< typename Element >
      uint32_t crc32(Element const *buffer, size_t length = 1);
  };
}

#include "ScidacChecksum/ScidacChecksum.inlines"
#include "ScidacChecksum/ScidacChecksum_aggregate.template"
#include "ScidacChecksum/ScidacChecksum_blockAggregate.template"
#include "ScidacChecksum/ScidacChecksum_crc32_a.template"
#include "ScidacChecksum/ScidacChecksum_crc32_b.template"
#include "ScidacChecksum/ScidacChecksum_calculate.template"

#endif
