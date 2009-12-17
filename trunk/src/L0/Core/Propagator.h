#pragma once

#include <L0/Base/Base.h>
#include <L0/Core/Field.h>
// #include <L0/QCD/Spinor.h>

namespace Core
{
  /* a propagator is a tensorial stucture with 2 Dirac and 2 Colour indices,
     each for source and sink, respectively, for all lattice sites 
     (being sink locations for a constant source)
     maybe at some point we should include the source properties, in whatever format,
     as a member variable
  */
  template< size_t L, size_t T >
  class Propagator
  {
    static const size_t nDirac  = 2;
    static const size_t nColour = 2;

    size_t psize;

    size_t *d_references;

    Core::Field< std::complex< double >, L, T > **d_components;

    size_t *colour_strides;
    size_t *dirac_strides;

    public:

// #include "Propagator/Propagator.view"
// #include "Propagator/Propagator.const_view"

      Propagator(bool alloc=true);
      Propagator(Propagator const &other);
//       Propagator(view< L, T > const &view); // this constructor makes a deep copy
      ~Propagator();

#include "Propagator/Propagator.operators"

#include "Propagator/Propagator.iterator"
// #include "Propagator/Propagator.const_iterator"

      iterator<1> begin(Base::ColourIndex const idx, const size_t ColourID);
      iterator<1> end(Base::ColourIndex const idx, const size_t ColourID);
      iterator<1> begin(Base::DiracIndex  const idx, const size_t DiracID);
      iterator<1> end(Base::DiracIndex  const idx, const size_t DiracID);

//       const_iterator begin(Base::ColourIndex const idx, const size_t ColourID) const;
//       const_iterator end(Base::ColourIndex const idx, const size_t ColourID) const;// 
//       const_iterator begin(Base::DiracIndex const idx, const size_t DiracID) const;
//       const_iterator end(Base::DiracIndex const idx, const size_t DiracID) const;
      

      size_t const size() const;
      size_t const numDirac() const;
      size_t const numColour() const;

/* NOTE we should include something like Propagator*Gamma*Propagator
   and maybe some herm. conj. derived class
   type cast of a propagator with one dirac and one colour index to a Core::Field< SU3::Spinor >
   would simplify I/O for standard ILDG format
*/ 

    private:
      void destroy();
      void isolate();

  };

}



// #include "Propagator/Propagator.view.inlines"
// #include "Propagator/Propagator.const_view.inlines"

#include "Propagator/Propagator.inlines"

#include "Propagator/Propagator.iterator.inlines"
// #include "Propagator/Propagator.const_iterator.inlines"

#include "Propagator/Propagator_destroy.template"
#include "Propagator/Propagator_isolate.template"


