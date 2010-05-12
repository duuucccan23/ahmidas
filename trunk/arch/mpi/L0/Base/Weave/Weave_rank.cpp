#include "Weave.ih"

size_t Base::Weave::rank(size_t const *coords) const
{
  int mpi_coord[4];
  for (size_t idx = 0; idx < 4; ++idx)
    mpi_coord[idx] = static_cast< int >(coords[idx]);
  return d_grid.rank(mpi_coord);
}
