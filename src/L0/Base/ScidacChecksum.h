#ifndef GUARD_BASE_SCIDACCHECKSUM_H
#define GUARD_BASE_SCIDACCHECKSUM_H

#include <cstdint>
#include <cstdlib>

namespace Base
{
  class ScidacChecksum
  {
    uint32_t d_sum29;
    uint32_t d_sum31;

    uint32_t d_crcTable[256];

    public:
      ScidacChecksum();
      void clear();

      template< typename Element >
      void aggregate(Element *data, size_t elements, size_t rank = 0);

    private:
      void createTable();
  };
}

#include "ScidacChecksum/ScidacChecksum.inlines"
//#include "ScidacChecksum/ScidacChecksum_aggregate.template"

#endif
