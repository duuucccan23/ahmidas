#ifndef CORE_UNIFOLD_H
#define CORE_UNIFOLD_H

#include <L0/Core/Grid.h>

// This particular version of Unifold is tailored for an MPI setup.

namespace Core
{
  template< size_t L, size_t T >
  class Unifold
  {
    static MPI::Grid< L, T > s_grid;

    public:
      static size_t const s_nodes;
      static size_t const s_volume;
      static size_t const s_size[4];
      static size_t const s_rank;
  };
}

template< size_t L, size_t T >
MPI::Grid Core::Unifold< L, T > s_grid();

// CHECK MPI DOCUMENTATION FOR THE FOLLOWING CALLS
template< size_t L, size_t T >
size_t const Core::Unifold::s_nodes = MPI::Get_size(s_grid.d_grid);

template< size_t L, size_t T >
size_t const Core::Unifold< L, T >::s_volume = L * L * L * T / s_nodes;

template< size_t L, size_t T >
size_t const Core::Unifold< L, T >::s_size[0] = ;

template< size_t L, size_t T >
size_t const Core::Unifold< L, T >::s_rank = MPI::Get_rank(s_grid.d_grid);

#endif
