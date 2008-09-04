#ifndef COM_H
#define COM_H

#include <L0/Core/Grid.h>

// MPI Communication handling

namespace Core
{
  template< typename Element, size_t L, size_t T >
  class Com
  {
    Grid< L, T>   &d_grid;
    MPI::Datatype  d_surfaces[4];

    public:
      Com();
      Com(Com const &other);
      Com< Element, L, T > &operator=(Com< Element, L, T > const &other);
//      ~Com();

  };

}

#include "Com/Com_Com.template" 

#endif
