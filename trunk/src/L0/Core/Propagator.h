#pragma once

#include <L0/Base/Base.h>
#include <L0/Core/Field.h>
// #include <L0/QCD/Spinor.h>

namespace Core
{
  // NOTE nDirac is the  n u m b e r  of Dirac indices, not a particualar Dirac index
  // same is true for nColour , being the number of colour indices
  template< size_t L, size_t T, size_t nDirac, size_t nColour >
  class Propagator
  {
    size_t *d_references;

    Core::Field< std::complex< double >, L, T > *const *d_components;

    size_t **colour_strides;
    size_t **dirac_strides;

    public:
      Propagator();
      Propagator(Propagator const &other);
      Propagator(view< L, T, nDirac, nColour > const &view);
      ~Propagator();

#include "Propagator/Propagator.operators"

#include "Propagator/Propagator.view"
// #include "Propagator/Propagator.const_view"

// #include "Propagator/Propagator.iterator"
// #include "Propagator/Propagator.const_iterator"

//       iterator begin(Base::ColourIndex const idx, const size_t ColourID);
//       iterator end(Base::ColourIndex const idx, const size_t ColourID);
//       iterator begin(Base::DiracIndex  const idx, const size_t DiracID);
//       iterator end(Base::DiracIndex  const idx, const size_t DiracID);

//       const_iterator begin(Base::ColourIndex const idx, const size_t ColourID) const;
//       const_iterator end(Base::ColourIndex const idx, const size_t ColourID) const;// 
//       const_iterator begin(Base::DiracIndex const idx, const size_t DiracID) const;
//       const_iterator end(Base::DiracIndex const idx, const size_t DiracID) const;

      view< L, T, nDirac-1, nColour > &operator()(Base::DiracIndex const idx, const size_t DiracID);
      view< L, T, nDirac, nColour-1 > &operator()(Base::ColourIndex const idx, const size_t ColourID);
      
//       const_view< L, T, nDirac-1, nColour > const &operator()(Base::DiracIndex const idx, const size_t DiracID) const;      
//       const_view< L, T, nDirac, nColour-1 > const &operator()(Base::ColourIndex const idx, const size_t ColourID) const;
      

      size_t size() const;

    private:
      void destroy();
      void isolate();
  };
}

#include "Propagator/Propagator.inlines"

#include "Propagator/Propagator.view.inlines"
// #include "Propagator/Propagator.const_view.inlines"

// #include "Propagator/Propagator.iterator.inlines"
// #include "Propagator/Propagator.const_iterator.inlines"

#include "Propagator/Propagator_destroy.template"
#include "Propagator/Propagator_isolate.template"
