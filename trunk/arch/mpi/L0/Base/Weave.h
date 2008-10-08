#ifndef GUARD_BASE_WEAVE_H
#define GUARD_BASE_WEAVE_H

#include <L0/Base/Grid.h>
#include <mpi/mpi.h>
// This particular version of Weave is tailored for an MPI setup.

namespace Base
{
  template< size_t L, size_t T >
  class Weave
  {
    static MPI::Grid<L,T> s_grid;

    public:
      static size_t const s_nodes;
      static size_t const s_volume;
      static size_t const s_size[4];
      static size_t const s_rank;
  };
}

template< size_t L, size_t T >
MPI::Grid Base::Weave< L, T > s_grid();

// CHECK MPI DOCUMENTATION FOR THE FOLLOWING CALLS
template< size_t L, size_t T >
size_t const Base::Weave::s_nodes = MPI::Get_size(s_grid.d_grid);

template< size_t L, size_t T >
size_t const Base::Weave< L, T >::s_volume = L * L * L * T / s_nodes;

template< size_t L, size_t T >
size_t const Base::Weave< L, T >::s_size[0] = ;

template< size_t L, size_t T >
size_t const Base::Weave< L, T >::s_rank = MPI::Get_rank(s_grid.d_grid);

#endif
