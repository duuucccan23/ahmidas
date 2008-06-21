#ifndef GUARD_CONTRACTION_MESON_H
#define GUARD_CONTRACTION_MESON_H

#include <Core/Field/Field.h>
#include <Dirac/Gamma/Gamma.h>
#include <QCD/Spinor/Spinor.h>

namespace Contraction
{
  template< size_t L, size_t T, size_t IndexSource, size_t IndexSink >
  class MesonLightLight
  {
    std::complex< double > *d_correlator;
    Core::Field< std::complex< double >, L, T > d_field;

    public:
      MesonLightLight(Core::Grid &grid, size_t sourceSlice);
      void contract(QCD::Spinor const &source, QCD::Spinor const &sink);

    private:
      void gatherCorrelator();
  };
}

#include "Meson.inlines"
#include "Meson_a.template"
#include "Meson_b.template"

#endif
