#pragma once

#include <vector>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Matrix.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>
#include <L0/Core/Component.h>
#include <L1/Smear.h>

namespace Contract
{

  Core::BaryonCorrelator proton_twopoint(Core::Propagator const &u1, Core::Propagator const &u2, Core::Propagator const &d,
                                   Base::BaryonPropagatorProjector const projector);

  Core::BaryonCorrelator proton_twopoint_alternative(Core::Propagator const &u1, Core::Propagator const &u2, 
                                                     Core::Propagator const &d,
                                                     Base::BaryonPropagatorProjector const projector);

  // for arbitrary interpolating field
  // G1a and G1b are the indices of the Gamma structures at source
  // G2a and G2b are associated to the sink Gamma structures
  template< size_t G1a, size_t G1b, size_t G2a, size_t G2b >
  Core::BaryonCorrelator nucleon_twopoint(Core::Propagator const &DoubleFlavourPhysBaseA, 
                                          Core::Propagator const &DoubleFlavourPhysBaseB,
                                          Core::Propagator const &SingleFlavorPhysBase,
                                          std::complex< double > const&globalFactor,
                                          Dirac::Gamma< G1a > const &GammaSourceA, 
                                          Dirac::Gamma< G1b > const &GammaSourceB, 
                                          Dirac::Gamma< G2a > const &GammaSinkA, 
                                          Dirac::Gamma< G2b > const &GammaSinkB); 

  /* ---------------------------------------------------------------------------------------------- */

  // keep full source and sink Dirac index dependence of the proton


  void create_sequential_source_proton_d(Core:: Propagator * const seqSrc,
                                         Core::Propagator const &u1, Core::Propagator const &u2,
                                         size_t const t_snk);

  void create_sequential_source_proton_u(Core::Propagator * const seqSrc,
                                         Core::Propagator const &d, Core::Propagator const &u,
                                         size_t const t_snk);


  /* ---------------------------------------------------------------------------------------------- */

  // standard versions, projector inserted at sink, trace is taken


  void create_sequential_source_proton_d(Core::Propagator &seqSrc,
                                         Core::Propagator const &u1, Core::Propagator const &u2,
                                         size_t const t_snk, Base::BaryonPropagatorProjector const projector);

  void create_sequential_source_proton_u(Core::Propagator &seqSrc,
                                         Core::Propagator const &d, Core::Propagator const &u,
                                         size_t const t_snk, Base::BaryonPropagatorProjector const projector);


  /* ---------------------------------------------------------------------------------------------- */

  // versions with smearing, expects (smeared) gauge field and source- and sink-smeared propagators.
  // the sequential source is smeared in order to have a fully smeared proton at the sink.


  void create_sequential_source_proton_d(Core:: Propagator &seqSrc,
                                        Core::Propagator const &u1, Core::Propagator const &u2,
                                        Core::Field < QCD::Gauge > &gauge_field,
                                        Smear::fermionFieldSmearing const smearing,
                                        // smearing: iterations and parameter
                                        size_t const nSmear, double const pSmear,
                                        size_t const t_snk, Base::BaryonPropagatorProjector const proj);

  void create_sequential_source_proton_u(Core:: Propagator &seqSrc,
                                         Core::Propagator const &d,
                                         Core::Propagator const &u,
                                         Core::Field < QCD::Gauge > &gauge_field,
                                         Smear::fermionFieldSmearing const smearing,
                                         // smearing: iterations and parameter
                                         size_t const nSmear, double const pSmear,
                                         size_t const t_snk, Base::BaryonPropagatorProjector const proj);


  // accordingly with free proton Dirac indices:


  void create_sequential_source_proton_d(Core:: Propagator * const seqSrc,
                                        Core::Propagator const &u1, Core::Propagator const &u2,
                                        Core::Field < QCD::Gauge > &gauge_field,
                                        Smear::fermionFieldSmearing const smearing,
                                        // smearing: iterations and parameter
                                        size_t const nSmear, double const pSmear,
                                        size_t const t_snk);

  void create_sequential_source_proton_u(Core:: Propagator * const seqSrc,
                                         Core::Propagator const &d,
                                         Core::Propagator const &u,
                                         Core::Field < QCD::Gauge > &gauge_field,
                                         Smear::fermionFieldSmearing const smearing,
                                         // smearing: iterations and parameter
                                         size_t const nSmear, double const pSmear,
                                         size_t const t_snk);


  /* ---------------------------------------------------------------------------------------------- */

  // yet another version that creates the sequential source at the operator insertion timeslice, such that
  // in the contraction, the sink timeslice is completely free

  void create_sequential_source_proton_fixed_insertion_timeslice(Core:: Propagator *seqSrc_u,
                                                                 Core:: Propagator *seqSrc_d,
                                                                 Core::Propagator const &u,
                                                                 Core::Propagator const &d,
                                                                 size_t const t_op, Base::Operator op);



  void create_sequential_source_proton_fixed_insertion_timeslice_derivative_operator(Core::Propagator *seqSrc_u,
										   Core::Propagator *seqSrc_d,
										   Core::Field< QCD::Gauge > &gauge_field,
										   Core::Propagator const &u,
										   Core::Propagator const &d,
										   size_t const t_op, Base::Operator op);


  /* ---------------------------------------------------------------------------------------------- */



  std::vector< Core::BaryonCorrelator > proton_threepoint_sequential(
    Core:: Propagator const &bw_prop_u, Core::Propagator const &fw_prop_u,
    Core:: Propagator const &bw_prop_d, Core::Propagator const &fw_prop_d,
    Core::Field< QCD::Gauge > * const gauge_field,
    std::vector< Base::Operator > ops);


  // does a set of contractions depending on the identifier
  /* identifier  --  operators
   * "local"         all Base::HermitianBilinearOperators in the order specified in Base.h
   * "1D"            gamma_mu D_nu (mu,nu=1,...,4, nu faster running index)
   *  -> addCorrs:   gamma_5 gamma_mu D_nu (mu,nu=1,...,4, nu faster running index)
   * "noether"       conserved vector current
   *  -> addCorrs:   conserved axial current
   */
  std::vector< Core::BaryonCorrelator > proton_threepoint_sequential(
    Core:: Propagator const &bw_prop_u, Core::Propagator const &fw_prop_u,
    Core:: Propagator const &bw_prop_d, Core::Propagator const &fw_prop_d,
    Core::Field< QCD::Gauge > * const gauge_field,
    std::string identifier, std::vector< Core::BaryonCorrelator > * addCorrs = NULL);



  /* ---------------------------------------------------------------------------------------------- */


  // the proton three point with a stochastic (estimate of an) all-to-all propagator at the proton sink.
  // this is a mess concerning performance since it scales vith 4-Volume^2, but provides a good cross-check
  // the gauge field is sometimes not needed, in this case NULL can be passed;
  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic_naive(Core::Propagator const &u,
                                                          Core::Propagator const &d,
                                                          Core::Propagator const &u_local,
                                                          Core::Propagator const &d_local,
                                                          Core::StochasticPropagator<12> const &u_stoch_at_sink,
                                                          Core::StochasticPropagator<12> const &d_stoch_at_sink,
                                                          Core::StochasticSource<12> const &xi_u_at_sink,
                                                          Core::StochasticSource<12> const &xi_d_at_sink,
                                                          Core::Field < QCD::Gauge > const * const gauge_field,
                                                          std::vector< Base::Operator > const &ops,
                                                          size_t const t_src, size_t const t_snk);

  // does a set of contractions depending on the identifier
  /* identifier  --  operators
   * "local"         all Base::HermitianBilinearOperators in the order specified in Base.h
   * "1D"            gamma_mu D_nu (mu,nu=1,...,4, nu faster running index)
   *  -> addCorrs:   gamma_5 gamma_mu D_nu (mu,nu=1,...,4, nu faster running index)
   * "noether"       conserved vector current
   *  -> addCorrs:   conserved axial current
   */
  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic(
                                                Core::Propagator const &S_u,
                                                Core::Propagator const &S_d,
                                                Core::Propagator const &S_u_local,
                                                Core::Propagator const &S_d_local,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi_u,
                                                Core::StochasticSource <12> const &xi_d,
                                                Core::Field<QCD::Gauge> const * const gauge_field,
                                                size_t const t_source, size_t const t_sink,
                                                std::string const &operator_flag,
                                                Base::BaryonPropagatorProjector const my_projector,
                                                std::vector< Core::BaryonCorrelator > &addCorrs);

  // same as above but without the projector (indices kept free)
  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic(
                                                Core::Propagator const &S_u,
                                                Core::Propagator const &S_d,
                                                Core::Propagator const &S_u_local,
                                                Core::Propagator const &S_d_local,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi_u,
                                                Core::StochasticSource <12> const &xi_d,
                                                Core::Field<QCD::Gauge> const * const gauge_field,
                                                size_t const t_source, size_t const t_sink,
                                                std::string const &operator_flag,
                                                std::vector< Core::BaryonCorrelator > &addCorrs);

  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic_non_local(Core::Propagator const &u,
                                                Core::Propagator const &d,
                                                Core::Propagator const &u_local,
                                                Core::Propagator const &d_local,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi_u,
                                                Core::StochasticSource <12> const &xi_d,
                                                Core::Field<QCD::Gauge> const * const gauge_field,
                                                size_t const t_source, size_t const t_sink,
                                                /* this allows for more than one operator */
                                                std::vector< Base::Operator > const &my_operators,
                                                std::vector< Base::BaryonPropagatorProjector > const &my_projectors_u,
                                                std::vector< Base::BaryonPropagatorProjector > const &my_projectors_d);
}

#include "Baryon/Baryon_nucleon_twopoint.template"

