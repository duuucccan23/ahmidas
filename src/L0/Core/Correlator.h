#pragma once

#include <vector>

#include <L0/Base/Weave.h>
#include <L0/Dirac/Matrix.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Tensor.h>

namespace Core
{
  template< typename Datatype >
  class Correlator;

  template< typename Datatype >
  std::ostream &operator<<(std::ostream &out, Correlator< Datatype > const &c);


  template< typename Datatype >
  class Correlator
  {

    friend class Weave;

    Base::Weave *d_weave;
    size_t *d_references;
    Field < Datatype > *d_data;
    Datatype *d_sumTimeslice;

    size_t d_T;
    size_t d_L;
    size_t d_offset;

    // for compatibility with parallel code
    Datatype *d_sumTimeslice_global;

	static Field< int * > *s_xRelative;

	  public:

	Correlator(Field < Datatype > *d_data);
	Correlator(Correlator const &other);

	~Correlator();

	Correlator &operator=(Correlator const &rhs);

	Datatype &operator[](size_t const idx);
	Datatype const &operator[](size_t const idx) const;


	void operator*=(double const factor);
	void operator*=(std::complex< double > const &factor);
	void operator*=(Base::BaryonPropagatorProjector const projector);
	void operator +=(Correlator< Datatype > const &other);
	void operator *=(Correlator< Datatype > const &other);

	void deleteField();

	// performs a summation over the volume of each individual timeslice (zero momentum projection)
	void sumOverSpatialVolume();
	// performs a summation over the volume of each individual timeslice including non-zero momentum projection
	std::vector< Correlator< Datatype > > momentumProjection(std::vector< int * > const &momenta) const;
	void prepareMomentumProjection(int const * const position_offset);

	void rotateToPhysicalBasis();

	// set offset (for conventional reasons one would often like to shift the source timeslice to 0)
	void setOffset(size_t timeslice);

	std::complex <double> getTrSum(size_t const timeslice) const;

	bool isRoot() const;

  // prints the trace and the momentum
  void printWithMomentum(std::ostream &out, int const * const momentum, std::string const& prefix="") const;
  // prints the full correlator and the momentum
  void printWithMomentum_full(std::ostream &out, int const * const momentum, std::string const& prefix="") const;

	size_t T() const;
	size_t L() const;
	size_t size() const;

	friend std::ostream &operator<< < Datatype >(std::ostream &out, Correlator< Datatype > const &c);

	private:
	 void destroy();
	 void isolate();
	 void isolate_action();
  };


#include "Correlator/Correlator.inlines"

#include "Correlator/Correlator_Correlator_a.template"
#include "Correlator/Correlator_Correlator_b.template"
#include "Correlator/Correlator.operators.template"

#include "Correlator/Correlator_deleteField.template"
#include "Correlator/Correlator_destroy.template"
#include "Correlator/Correlator_isolate.template"
#include "Correlator/Correlator_momentumProjection.template"
#include "Correlator/Correlator_prepareMomentumProjection.template"
#include "Correlator/Correlator_setOffset.template"
#include "Correlator/Correlator_sumOverSpatialVolume.template"

  typedef Correlator< Dirac::Matrix > BaryonCorrelator;

}
