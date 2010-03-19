#pragma once

#include <vector>
#include <string>
#include <L0/Base/Base.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/Core/Correlator.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Tensor.h>

namespace Core
{
  /* a propagator is a tensorial stucture with 2 Dirac and 2 Colour indices,
     each for source and sink, respectively, for all lattice sites 
     (being sink locations for a constant source)
     maybe at some point we should include the source properties, in whatever format,
     as a member variable
  */
  template< size_t NComp >
  class StochasticSource;

  template< size_t NComp >
  class StochasticPropagator;

  class Propagator
  {

    size_t *d_references;

    static const size_t nDirac  = 2;
    static const size_t nColour = 2;

    static const size_t d_size = 144;

    protected:
    Core::Field< QCD::Tensor > *d_components;


    public:

      Propagator(size_t L, size_t T, bool alloc=true);
      Propagator(Propagator const &other);
      ~Propagator();


      //Propagator &operator=(Propagator const &rhs);

      // return entry at particular lattice site, corresponding to the propagator sink
      // for parallelization reasons this does not return a reference
      QCD::Tensor operator()(size_t const* sinkSite) const;

      template< size_t Index >
      void operator*=(Dirac::Gamma< Index > const &gamma);

      template< size_t Index >
      Propagator operator*(Dirac::Gamma< Index > const &gamma) const;

      // needed for meson contractions
      Core::Field< QCD::reducedTensor > *operator*(Propagator const &other) const;

      // needed for baryon contractions
      Core::Field< QCD::reducedTensor > *construct_baryon(Propagator const &no2, Propagator const &no3,
                                                          Base::BaryonInterpolatingField const ipol) const;

      /*
          Revert Propagator using gamma5 hermeticity trick:
          S (x,y) = gamma5 * S^dagger (y,x) * gamma5,
          where the dagger is meant in Dirac and colour space only.
          Note:
           - for point or similar source x is generally fixed
           - for twisted mass this also changes the flavour
      */
      Propagator &revert();

      Propagator &dagger();

      // average difference of two different propagators
      double diff(Propagator const& other) const;

      // just for testing
      void setToRandom();

      void changeBoundaryConditions_uniformToFixed(size_t timesliceSource, size_t timesliceBoundary);

      Propagator &shift(Base::SpaceTimeIndex const idx, Base::Direction const dir, int const times);

      void gaugeTransform_fixedSource(Field< SU3::Matrix > const &gaugeTrafo, size_t const * sourcePos);


#include "Propagator/Propagator.iterator"
#include "Propagator/Propagator.const_iterator"

      iterator begin();
      iterator end();
      const_iterator begin() const;
      const_iterator end() const;

      size_t const size() const;
      size_t const L() const;
      size_t const T() const;

      template< size_t Index >
      friend Propagator operator*(Dirac::Gamma< Index > const &gamma, Propagator const &p);

      friend std::ostream &operator<<(std::ostream &out, Propagator const &p);

    private:
      void destroy();
      void isolate();

  };

  template< size_t Index >
  Propagator operator*(Dirac::Gamma< Index > const &gamma, Propagator const &p);

  std::ostream &operator<<(std::ostream &out, Propagator const &p);

  template< size_t NComp >
  class StochasticSource : public Propagator
  {

    friend class StochasticPropagator< NComp >;

    public:

      StochasticSource< NComp > (size_t const L, size_t const T);
      StochasticSource< NComp > (size_t const L, size_t const T,
                                 Base::SourcePolarization const, Base::SourceColorState const);

      explicit StochasticSource< NComp > (Propagator const &base);
      StochasticSource< NComp > (StochasticSource< NComp > const &other);

      // ~StochasticSource< NComp > ();

      // void conjugate();

      Propagator operator*(StochasticPropagator< NComp > const &sPropagator) const;

      Propagator createStochasticPropagator_fixedSink(StochasticPropagator< NComp > const &, size_t const *) const;

      size_t const L() const;
      size_t const T() const;

  };

  template< size_t NComp >
  class StochasticPropagator : public Propagator
  {

    friend class StochasticSource< NComp >;

    public:

      StochasticPropagator< NComp > (size_t const L, size_t const T);

      StochasticPropagator< NComp > (Propagator const &other);

      StochasticPropagator< NComp > (StochasticPropagator< NComp > const &other);

      Field< QCD::reducedTensor > *operator*(StochasticPropagator< NComp > const &other) const;

      Propagator operator*(StochasticSource< NComp > const &sSource) const;

  };

}


#include "Propagator/Propagator.inlines"
#include "Propagator/Propagator.baryon.inlines"
#include "Propagator/StochasticPropagator.inlines"
#include "Propagator/StochasticSource.inlines"

#include "Propagator/Propagator.iterator.inlines"
#include "Propagator/Propagator.const_iterator.inlines"

#include "Propagator/Propagator.operators.inlines"

#include "Propagator/Propagator_destroy.template"
#include "Propagator/Propagator_isolate.template"




