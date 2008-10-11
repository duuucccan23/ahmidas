#ifndef GUARD_SOURCE_SLICE_H
#define GUARD_SOURCE_SLICE_H

namespace Source
{
  template< size_t L, size_t T >
  class Slice
  {
    Field< QCD::Spinor, L, 1 >  d_source;

    public:
      Slice(Stochastic const &source, Base::ColourIndex, Base::DiracIndex);
  };
}

#endif
