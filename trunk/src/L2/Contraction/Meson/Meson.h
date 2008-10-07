#ifndef GUARD_CONTRACTION_MESON_H
#define GUARD_CONTRACTION_MESON_H

#include <Base/Grid.h>
#include <Core/Field.h>
#include <Dirac/Gamma.h>
#include <QCD/Spinor.h>

namespace Contraction
{
  template< size_t L, size_t T, size_t IndexSource, size_t IndexSink >
  class Meson
  {
    std::complex< double >                      *d_correlator;
    Core::Field< std::complex< double >, L, T >  d_field;

    public:
      Meson(Base::Grid &grid);
      void contract(Core::Field< QCD::Spinor, L, T > const &source, Core::Field< QCD::Spinor, L, T > QCD::Spinor const &sink);

    private:
      void gatherCorrelator();
  };
}

#include "Meson.inlines"
#include "Meson_a.template"

#endif
