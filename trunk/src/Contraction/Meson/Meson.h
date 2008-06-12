#ifndef GUARD_CONTRACTION_MESON_H
#define GUARD_CONTRACTION_MESON_H

#include "../../Core/Field/Field.h"
#include "../../Dirac/Gamma/Gamma.h"
#include "../../QCD/Spinor/Spinor.h"

namespace Contraction
{
  template< size_t IndexSource, size_t IndexSink >
  class MesonLightLight
  {
    Core::Grid &d_grid;

    std::complex< double > *d_correlator;
    std::complex< double > *d_buffer;

    public:
      MesonLightLight(Core::Grid &grid, size_t sourceSlice);
      void contract(QCD::Spinor const &source, QCD::Spinor const &sink);

    private:
      void gatherCorrelator();
  };
}

#include "MesonLightLight.inlines"

#endif
