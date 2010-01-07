#include "Weave.ih"

Base::Weave::Weave(Weave const &other)
: d_L(other.d_L), d_T(other.d_T), d_localVolume(other.d_localVolume), d_globalVolume(other.d_globalVolume)
{
  std::copy(other.d_surfaces, other.d_surfaces + 4, d_surfaces);
}
