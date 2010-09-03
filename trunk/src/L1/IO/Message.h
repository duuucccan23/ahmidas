#pragma once

#include <string>

namespace IO
{
  namespace Lime
  {
    class Message
    {
      size_t      d_message;
      size_t      d_record;
      std::string d_type;
      size_t      d_datalen;
      size_t      d_padlen;
      bool        d_MBflag;
      bool        d_MEflag;
      std::string d_data;
      bool        d_isXML;

      public:
        Message();

    }
  }
}