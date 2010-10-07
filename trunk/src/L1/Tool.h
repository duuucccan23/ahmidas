#pragma once

#include <vector>
#include <string>

#include <L0/Base/Base.h>
#include <L0/Core/Field.h>
#include <L0/Core/Correlator.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/SU3/Matrix.h>

namespace Tool
{
  void reunitarize(Core::Field< QCD::Gauge > *field);

  void reunitarize(Core::Field< SU3::Matrix > *field);

  template< typename Element >
  void randomize(Core::Field< Element > *field);

//   void randomize(Core::Propagator *propagator);

  // print correlation functions in a quasi-standard ETMC format
  void printLightMesonCorrelator(std::vector< Core::Correlator< Dirac::Matrix > > const &correlator,
                                 std::string const filename = "output.txt");

  template< typename Element >
  void setToZero(Core::Field< Element > *field);

  template< typename Element >
  void setToIdentity(Core::Field< Element > *field);

  std::complex< double > tr(Core::Field< SU3::Matrix > const &field);
  double realtr(Core::Field< SU3::Matrix > const &field);
  Core::Field< std::complex< double > > localTrace(Core::Field< SU3::Matrix > const &field);
  Core::Field< double > localRealTrace(Core::Field< SU3::Matrix > const &field);

  double spatialPlaquette(Core::Field< QCD::Gauge > const &field);
  double spatialDownPlaquette(Core::Field< QCD::Gauge > const &field);
  double temporalPlaquette(Core::Field< QCD::Gauge > const &field);
  double temporalDownPlaquette(Core::Field< QCD::Gauge > const &field);
  void fixCoulombGauge(Core::Field< QCD::Gauge > *field);

  SU3::Matrix killTrace(SU3::Matrix const &target);

  //Please not that in the following, the left spinor field is complex conjugated automatically in this inner product!
  std::complex < double > innerProduct(Core::Field< QCD::Spinor > const &left,
                                       Core::Field< QCD::Spinor > const &right,
                                       const size_t tslice);

  Base::Operator convertIntToOperator(int const input);

}

#include "Tool/Tool_randomize.template"
#include "Tool/Tool_setToZero.template"
#include "Tool/Tool_setToIdentity.template"
