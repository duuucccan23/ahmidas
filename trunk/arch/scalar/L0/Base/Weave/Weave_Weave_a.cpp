#include "Weave.ih"

Base::Weave::Weave(size_t const L, size_t const T)
: d_L(L), d_T(T), d_localVolume(L * L * L * T), d_globalVolume(d_localVolume)
{
  d_surfaces[idx_X] = 1;
  d_surfaces[idx_Y] = L * d_surfaces[idx_X];
  d_surfaces[idx_Z] = L * d_surfaces[idx_Y];
  d_surfaces[idx_T] = L * d_surfaces[idx_Z];
}
