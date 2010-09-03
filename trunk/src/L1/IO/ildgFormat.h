#pragma once

#include <string>
#include <L1/IO/Lime/Message.h>

namespace IO
{
  namespace Lime
  {
    class ildgFormat : public Message
    {
      static std::string s_start = "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"        xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">";
      static std::string s_end = "</ildgFormat>";
      size_t      d_nx;
      size_t      d_ny;
      size_t      d_nz;
      size_t      d_nt;
      size_t      d_precision;
      std::string d_field;
      std::string d_version;

      public:
        ildgFormat();

    }
  }
}