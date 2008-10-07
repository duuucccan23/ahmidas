#ifndef GUARD_CORE_PROPAGATOR_H
#define GUARD_CORE_PROPAGATOR_H

#include <L0/Core/Core.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>

namespace Core
{
  template< size_t L, size_t T >
  class Propagator
  {
    Core::Field< Spinor, L, T > d_components[12];

    public:
      Propagator(Propagator const &other);
      Propagator(Core::Field< Spinor, L, T > *data);

#include "Propagator/Propagator.iterator"
#include "Propagator/Propagator.const_iterator"

      iterator_full begin();
      iterator_full end();

      const_iterator_full begin() const;
      const_iterator_full end() const;

      iterator_colour begin(Core::ColourIndex const idx);
      iterator_colour end(Core::ColourIndex const idx);

      const_iterator_colour begin(Core::ColourIndex const idx) const;
      const_iterator_colour end(Core::ColourIndex const idx) const;

      iterator_dirac begin(Core::DiracIndex const idx);
      iterator_dirac end(Core::DiracIndex const idx);

      const_iterator_dirac begin(Core::DiracIndex const idx) const;
      const_iterator_dirac end(Core::DiracIndex const idx) const;

#include "Propagator/Propagator.operators"
  };
}

#include "Propagator/Propagator.inlines"
#include "Propagator/Propagator.iterator.inlines"
#include "Propagator/Propagator.const_iterator.inlines"

#endif
