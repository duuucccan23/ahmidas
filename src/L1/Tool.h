#pragma once

#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
//#include <L0/Core/Propagator.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>
#include <L1/Path.h>

namespace Tool
{
  template< size_t L, size_t T >
  void reunitarize(Core::Field< QCD::Gauge, L, T > *field);

  template< size_t L, size_t T >
  void reunitarize(Core::Field< SU3::Matrix, L, T > *field);

  template< typename Element, size_t L, size_t T >
  void randomize(Core::Field< Element, L, T > *field);

//   template< size_t L, size_t T >
//   void randomize(Core::Propagator< L, T > *propagator);

  template< typename Element, size_t L, size_t T >
  void setToZero(Core::Field< Element, L, T > *field);

  template< typename Element, size_t L, size_t T >
  void setToIdentity(Core::Field< Element, L, T > *field);

  template< typename Element, size_t L, size_t T >
  std::complex< double > tr(Core::Field< Element, L, T > const &field);

  template< size_t L, size_t T >
  Core::Field< SU3::Matrix, L, T > localTrace(Core::Field< SU3::Matrix, L, T > const &field);

  template< size_t L, size_t T >
  double spatialPlaquette(Core::Field< QCD::Gauge, L, T > &field);

  template< size_t L, size_t T >
  double temporalPlaquette(Core::Field< QCD::Gauge, L, T > &field);

  template< size_t L, size_t T >
  void fixCoulombGauge(Core::Field< QCD::Gauge, L, T > *field);

  SU3::Matrix killTrace(SU3::Matrix const &target);
}

#include "Tool/Tool.inlines"

#include "Tool/Tool_reunitarize.template"
#include "Tool/Tool_randomize.template"
#include "Tool/Tool_setToZero.template"
#include "Tool/Tool_localTrace.template"
#include "Tool/Tool_spatialPlaquette.template"
#include "Tool/Tool_temporalPlaquette.template"
#include "Tool/Tool_fixCoulombGauge.template"
