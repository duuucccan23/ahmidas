#pragma once

#include <L1/IO/Param.h>

namespace IO
{
  class scidacChecksum : public Param
  {
    static char *s_start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacChecksum>";
    static char *s_end = "</scidacChecksum>";
    public :
      size_t       suma;
      size_t       sumb;
      std::string *d_version;

      virtual void parse();
      virtual void generate() const;
      char const *data() const;
  }
}

