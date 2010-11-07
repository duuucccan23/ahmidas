#pragma once

#include <vector>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/Core/Correlator.h>
#include <L0/QCD/Gauge.h>
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

    static const size_t nDirac  = 2;
    static const size_t nColour = 2;
    static const size_t d_size = 144;

    protected:
      Core::Field< QCD::Tensor > *d_components;

      size_t *d_references;

    public:

      Propagator(size_t L, size_t T, bool alloc=true);
      Propagator(Propagator const &other);
      ~Propagator();


      // return entry at particular lattice site, corresponding to the propagator sink
      // for parallelization reasons this does not return a reference
      QCD::Tensor operator()(size_t const* sinkSite) const;

      // unsafe access of Field elements
      QCD::Tensor       &operator[](size_t const localIndex);
      QCD::Tensor const &operator[](size_t const localIndex) const;

      Propagator &operator=(Propagator const &rhs);

      void operator*=(std::complex< double > const &factor);

      void operator+=(Propagator const &other);
      void operator-=(Propagator const &other);

      template< size_t Index >
      void operator*=(Dirac::Gamma< Index > const &gamma);

      template< size_t Index >
      Propagator operator*(Dirac::Gamma< Index > const &gamma) const;

      template < size_t Index >
      Propagator &rightMultiply(Dirac::Gamma< Index > const& gamma);

      // multiplication by operator (this is multiplication from the right)
      void rightMultiplyOperator(Base::Operator const O); // only for local operators
      void leftMultiplyOperator(Base::Operator const O); // only for local operators
      Core::Field< Dirac::Matrix > *contractWithOperatorInsertion(Base::Operator const O, Field< QCD::Gauge > * const gauge_field, Propagator const &fromRight);

      void rightMultiply(Propagator const & other);
      void leftMultiply(Propagator const & other);

      // Dirac operators:
      // fixed boundary conditions (light doublet)
      // the boundary conditions (if no theta-parameter is given) are periodic in space and antiperiodic in time
      Propagator applyDiracOperator(Field< QCD::Gauge > const &gauge_field, double const kappa, double const mu,
                                    double const thetaT = 1.0, double const thetaX = 0.0,
                                    double const thetaY = 0.0, double const thetaZ = 0.0) const;
      // uniform boundary conditions (heavy doublet)
      Propagator applyDiracOperator(Field< QCD::Gauge > const &gauge_field, double const kappa, double const m0, double const mu_sigma, double const mu_delta, Propagator const &secondProp, bool sign) const;
      // fixed boundary conditions (light doublet)
      Propagator applyDiracOperator(Field< QCD::Gauge > const &gauge_field, double const kappa, double const mu, size_t const t_boundary) const;

      // this function reconstructs the doublet propagators from  the (Q^+Q^-)^-1 output of the tmLQCD multi-mas-solver
      // the boundary conditions (if no theta-parameter is given) are periodic in space and antiperiodic in time
      void reconstruct_doublet(Propagator &propPlus, Propagator &propMinus,
                              Field< QCD::Gauge > const &gauge_field,
                              double const kappa, double const mu, double const thetaT = 1.0,
                              double const thetaX = 0.0, double const thetaY = 0.0, double const thetaZ = 0.0) const;

      Propagator &select_timeslice(size_t const timeslice);

      Propagator &smearJacobi(double const alpha, size_t const iterations, Field< QCD::Gauge > &gauge_field);
      Propagator &smearJacobi(double const alpha, size_t const iterations, Field< QCD::Gauge > &gauge_field, size_t const timeslice);

      // needed for meson contractions
      Core::Field< Dirac::Matrix > *operator*(Propagator const &other) const;

      // this is a simple contraction: A be a Tensor of *this, B be a Tensor of other,
      // both at the same lattice site.
      // this routine does A.leftMultiply(B) and takes the color trace of the result.
      Core::Field< Dirac::Matrix > *contract(Propagator const &other) const;

      // needed for baryon contractions (twopoint)
      Core::Field< Dirac::Matrix > *construct_baryon(Propagator const &no2, Propagator const &no3,
                                                          Base::BaryonInterpolatingField const ipol) const;

      // needed for baryon contractions (threepoint)
      std::vector< Core::Field< Dirac::Matrix > * >
        construct_baryon_with_operator_insertion(Propagator const &no2, Propagator const &no3,
                                                 StochasticPropagator< 12 > const &phi_no1,
                                                 StochasticPropagator< 12 > const &phi_no2,
                                                 StochasticPropagator< 12 > const &phi_no3,
                                                 StochasticSource< 12 > const &xi,
                                                 Base::BaryonInterpolatingField const ipol,
                                                 std::vector< Base::Operator > const &ops,
                                                 size_t const t_src, size_t const t_snk) const;


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
      Propagator &conjugate();
      Propagator &transpose();


      Field< std::complex<double> > trace() const;


      // average difference of two different propagators
      double diff(Propagator const& other) const;

      // just for testing
      void setToRandom();

      void changeBoundaryConditions_uniformToFixed(size_t timesliceSource, size_t timesliceBoundary);

      // sign = true  <=> positive tau_3 entry
      // sign = false <=> negative tau_3 entry
      void rotateToPhysicalBasis(bool const sign = true);


      Propagator &shift(Base::SpaceTimeIndex const idx, Base::Direction const dir, size_t const times=1);

      void gaugeTransform_fixedSource(Field< SU3::Matrix > const &gaugeTrafo, size_t const * sourcePos);


#include "Propagator/Propagator.iterator"
#include "Propagator/Propagator.const_iterator"

      iterator begin();
      iterator end();
      const_iterator begin() const;
      const_iterator end() const;

      double norm() const;
      double normq() const;
      size_t const size() const;
      size_t const L() const;
      size_t const T() const;

      template< size_t Index >
      friend Propagator operator*(Dirac::Gamma< Index > const &gamma, Propagator const &p);

      friend std::ostream &operator<<(std::ostream &out, Propagator const &p);

      void isolate();

    protected:
      void destroy();

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
                                 Base::SourcePolarization const pol, Base::SourceColorState const col,
                                 size_t const timeslice,
                                 Field< uint64_t > const &seeds,
                                 Base::SourceStochasticTypeFlag const type);

      explicit StochasticSource< NComp > (Propagator const &base);
      StochasticSource< NComp > (StochasticSource< NComp > const &other);

      // ~StochasticSource< NComp > ();

      //void conjugate();

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

      Field< Dirac::Matrix > *operator*(StochasticPropagator< NComp > const &other) const;

      StochasticPropagator< NComp > &select_timeslice(size_t const timeslice);

      Propagator operator*(StochasticSource< NComp > const &sSource) const;

      StochasticPropagator< NComp > &revert();
      StochasticPropagator< NComp > &dagger();
      StochasticPropagator< NComp > &conjugate();
      StochasticPropagator< NComp > &transpose();

      void isolate();

      void rotateToPhysicalBasis(bool const sign);
  };

  #include "Propagator/Propagator.inlines"
  #include "Propagator/StochasticPropagator.inlines"
  #include "Propagator/Propagator.iterator.inlines"
  #include "Propagator/Propagator.const_iterator.inlines"
  #include "Propagator/Propagator.operators.inlines"
}
#include "Propagator/StochasticSource.template"

