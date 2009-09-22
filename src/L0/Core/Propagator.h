#pragma once

#include <L0/Base/Base.h>
#include <L0/Core/Field.h>
// #include <L0/QCD/Spinor.h>

namespace Core
{
  // NOTE Dirac is the  n u m b e r  of Dirac indices, not a particualar Dirac index
  // same is true for Colour 
  template< size_t L, size_t T, size_t Dirac, size_t Colour >
  class Propagator
  {
    size_t *d_references;

    Core::Field< std::complex< double >, L, T > *d_components;

    size_t *colour_strides;
    size_t *dirac_strides;

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

      iterator_colour begin(Base::ColourIndex const idx, const size_t ColourID);
      iterator_colour end(Base::ColourIndex const idx, const size_t ColourID);

      const_iterator_colour begin(Base::ColourIndex const idx, const size_t ColourID) const;
      const_iterator_colour end(Base::ColourIndex const idx, const size_t ColourID) const;

      iterator_dirac begin(Base::DiracIndex const idx, const size_t DiracID);
      iterator_dirac end(Base::DiracIndex const idx, const size_t DiracID);

      const_iterator_dirac begin(Base::DiracIndex const idx, const size_t DiracID) const;
      const_iterator_dirac end(Base::DiracIndex const idx, const size_t DiracID) const;

      Propagator< L, T, Dirac-1, Colour > &getDirac(Base::DiracIndex const dirIdx, const size_t DiracID);
      Propagator< L, T, Dirac-1, Colour > const &getDirac(Base::DiracIndex const dirIdx, const size_t DiracID) const;

      Propagator< L, T, Dirac, Colour-1 > &getColour(Base::ColourIndex const dirIdx, const size_t ColourID);
      Propagator< L, T, Dirac, Colour-1 > const &getColour(Base::ColourIndex const dirIdx, const size_t ColourID) const;

      size_t size() const;

    private:
      void destroy();
      void isolate();
      void getStrides();
  };
}

#include "Propagator/Propagator.inlines"
#include "Propagator/Propagator.iterator.inlines"
#include "Propagator/Propagator.const_iterator.inlines"

#include "Propagator/Propagator_destroy.template"
#include "Propagator/Propagator_isolate.template"
