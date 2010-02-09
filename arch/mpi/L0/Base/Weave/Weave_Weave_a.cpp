#include "Weave.ih"

Base::Weave::Weave(size_t L, size_t T)
: d_L(L), d_T(T), d_grid(Base::Grid(L, T))
{
  d_globalVolume = d_L*d_L*d_L*d_T;
  d_localVolume = d_grid.localVolume();
  for (size_t idx = 0; idx < 4; ++idx)
  {
    d_surfaces[idx] = d_grid.surface(idx); //Surface size in direction
    d_localSize[idx] = d_grid.size(idx); //Local size in direction
  }
}
