#pragma once

#include <vector>
#include <string>
#include <L0/Base/Base.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Tensor.h>
#include <L0/Tool/IO.h>


namespace Core
{
  /* a propagator is a tensorial stucture with 2 Dirac and 2 Colour indices,
     each for source and sink, respectively, for all lattice sites 
     (being sink locations for a constant source)
     maybe at some point we should include the source properties, in whatever format,
     as a member variable
  */
  class Propagator
  {

    size_t *d_references;

    Core::Field< QCD::Tensor > *d_components;

    enum PropagatorStride
    {
      ColourStrideSink   =  1,
      ColourStrideSource = 12,
      DiracStrideSink    =  3,
      DiracStrideSource  = 36
    };

    static const size_t nDirac  = 2;
    static const size_t nColour = 2;
    size_t *colour_strides;
    size_t *dirac_strides;

    static const size_t d_size = 144;


    public:

      Propagator(size_t L, size_t T, bool alloc=true);
      Propagator(Propagator const &other);
      ~Propagator();

      bool load(std::vector< std::string > const filenames, std::string const format);


      template< size_t Index >
      void operator*=(Dirac::Gamma< Index > const &);


#include "Propagator/Propagator.iterator"
// #include "Propagator/Propagator.const_iterator"

      iterator begin(size_t const timeslice);
      iterator end(size_t  const timeslice);

      size_t const size() const;

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


#include "Propagator/Propagator.inlines"
#include "Propagator/Propagator.operators.inlines"

#include "Propagator/Propagator.iterator.inlines"
// #include "Propagator/Propagator.const_iterator.inlines"

#include "Propagator/Propagator_destroy.template"
#include "Propagator/Propagator_isolate.template"


