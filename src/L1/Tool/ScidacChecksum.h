#pragma once

#include <L0/Core/Field.h>

#include <cstdlib>
#include <stdint.h>

namespace Tool
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
      ScidacChecksum(uint64_t const &init);
      ScidacChecksum(ScidacChecksum const &other);
      void clear();

      uint64_t checksum() const;
      uint32_t lower() const;
      uint32_t upper() const;

      template< typename Element >
      size_t aggregate(Element const *data, size_t elements, size_t rank);

      template< typename Element >
      size_t aggregate(Element const &data, size_t rank);

      template< typename Element >
      size_t blockAggregate(Element const *data, size_t blockSize, size_t blocks, size_t rank);

      template< typename Element >
      size_t blockAggregate(Element const *data, Element const *finish, size_t blockSize, size_t rank);

      template< typename Element >
      size_t blockAggregate(Element const &data, size_t blockSize, size_t blocks, size_t rank);

      template< typename Element >
      void calculate(Core::Field< Element > const &field);

    public:
      void createTable();

      template< typename Element >
      uint32_t crc32(Element const &buffer, uint32_t crc);

      template< typename Element >
      uint32_t crc32(Element const *buffer, size_t length, uint32_t crc);
  };
}

#include "ScidacChecksum/ScidacChecksum.inlines"
#include "ScidacChecksum/ScidacChecksum_aggregate_a.template"
#include "ScidacChecksum/ScidacChecksum_aggregate_b.template"
#include "ScidacChecksum/ScidacChecksum_blockAggregate_a.template"
#include "ScidacChecksum/ScidacChecksum_blockAggregate_b.template"
#include "ScidacChecksum/ScidacChecksum_blockAggregate_c.template"
#include "ScidacChecksum/ScidacChecksum_crc32_a.template"
#include "ScidacChecksum/ScidacChecksum_crc32_b.template"
#include "ScidacChecksum/ScidacChecksum_calculate.template"
