#include "Weave.ih"

Base::Weave::Weave()
: d_grid(Base::Grid())
{
  d_localVolume = d_grid.localVolume();
  for (size_t idx = 0; idx < 4; ++idx)
  {
    d_surfaces[idx] = d_grid.surface(idx); //Surface size in direction
    d_localSize[idx] = d_grid.size(idx); //Local size in direction
  }
}
