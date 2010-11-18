#include "Weave.ih"

Base::Weave::Weave(size_t L, size_t T)
: d_L(L), d_T(T)
{
  if(s_grid_ref == 0)
  {
    s_grid = reinterpret_cast< void * >(new Base::Grid(L,T));
    d_grid = reinterpret_cast< Base::Grid *>(s_grid);
  }
  else
  {
    assert(s_grid != NULL);
    d_grid = reinterpret_cast< Base::Grid * >(s_grid);
    assert(d_grid->L() == d_L && d_grid->T() == d_T);
  }
  s_grid_ref++;
  d_globalVolume = d_L * d_L * d_L * d_T;
  d_localVolume = d_grid->localVolume();
  for (size_t idx = 0; idx < 4; ++idx)
  {
    d_surfaces[idx] = d_grid->dimSize(idx); //Surface size in direction
    d_localSize[idx] = d_grid->size(idx); //Local size in direction
  }
}
