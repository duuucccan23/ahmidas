#include "Weave.ih"

Base::Weave &Base::Weave::operator=(Base::Weave const &other)
{
  if (&other == this)
    return *this;

  std::copy(other.d_surfaces, other.d_surfaces + 4, d_surfaces);
  d_L = other.d_L;
  d_T = other.d_T;
  d_localVolume = other.d_localVolume;

  return *this;
}
