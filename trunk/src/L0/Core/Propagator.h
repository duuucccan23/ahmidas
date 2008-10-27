#ifndef GUARD_CORE_PROPAGATOR_H
#define GUARD_CORE_PROPAGATOR_H

#include <L0/Base/Base.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>

namespace Core
{
  template< size_t L, size_t T >
  class Propagator
  {
    size_t                            *d_references;
    Core::Field< QCD::Spinor, L, T > **d_components;

    public:
      Propagator();
      Propagator(Propagator const &other);
      ~Propagator();

#include "Propagator/Propagator.operators"

#include "Propagator/Propagator.iterator"
#include "Propagator/Propagator.const_iterator"

      iterator_full begin();
      iterator_full end();

      const_iterator_full begin() const;
      const_iterator_full end() const;

      iterator_colour begin(Base::ColourIndex const idx);
      iterator_colour end(Base::ColourIndex const idx);

      const_iterator_colour begin(Base::ColourIndex const idx) const;
      const_iterator_colour end(Base::ColourIndex const idx) const;

      iterator_dirac begin(Base::DiracIndex const idx);
      iterator_dirac end(Base::DiracIndex const idx);

      const_iterator_dirac begin(Base::DiracIndex const idx) const;
      const_iterator_dirac end(Base::DiracIndex const idx) const;

      size_t size() const;

    private:
      void destroy();
      void isolate();
  };
}

#include "Propagator/Propagator.inlines"
#include "Propagator/Propagator.iterator.inlines"
#include "Propagator/Propagator.const_iterator.inlines"

#include "Propagator/Propagator_destroy.template"
#include "Propagator/Propagator_isolate.template"

#endif
