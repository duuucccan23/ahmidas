#ifndef GUARD_BASE_IO_H
#define GUARD_BASE_IO_H

#include <string>
#include <L0/Tool/IO/Lime/Reader.h>
#include <L0/Tool/IO/Lime/Writer.h>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>

// Below is necessary to prevent circular including when we make IO functions friends of field.

namespace Base
{
  namespace IO
  {
    struct ILDGinfo
    {
      std::string version;
      std::string field;
      std::string precision;
      size_t      dims[4];

      ILDGinfo(Lime::Reader &reader);
    };

    template< typename Element, size_t L, size_t T >
    Core::Field< Element, L, T > loadILDG(std::string const &filename);

    template< typename Element, size_t L, size_t T >
    Core::Field< Element, L, T > loadScidac(std::string const &filename);

    template< typename Element, size_t L, size_t T >
    void saveILDG(Core::Field< Element, L, T > const &field, std::string const &filename);

    template< typename Element, size_t L, size_t T >
    void saveScidac(Core::Field< Element, L, T > const &field, std::string const &filename);
  }
}
#include "IO/loadILDG.template"
#include "IO/loadScidac.template"
#include "IO/saveILDG.template"
#include "IO/saveScidac.template"

#endif
