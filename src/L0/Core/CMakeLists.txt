SET(Sources_Correlator Correlator/Correlator_sumOverSpatialVolume.cpp Correlator/Correlator_projection.cpp
                       Correlator/Correlator_cout_operator_lshift.cpp Correlator/Correlator_isolate.cpp
                       Correlator/Correlator_destroy.cpp)
SET(Sources_Propagator Propagator/StochasticPropagator_StochasticPropagator.cpp
                       Propagator/StochasticSource.cpp
                       Propagator/Propagator_changeBoundaryConditions.cpp
                       Propagator/Propagator_construct_baryon.cpp
                       Propagator/Propagator_construct_baryon_with_operator_insertion.cpp
                       Propagator/Propagator_contract_a.cpp
                       Propagator/Propagator_diff.cpp
                       Propagator/Propagator_gaugeTransformation.cpp
                       Propagator/Propagator_isolate.cpp
                       Propagator/Propagator_Propagator.cpp
                       Propagator/Propagator_revert.cpp
                       Propagator/Propagator_shift.cpp
                       Propagator/Propagator_smearJacobi.cpp
                       Propagator/Propagator.operators.cpp)
ADD_LIBRARY(Core ${Sources_Correlator} ${Sources_Propagator})
TARGET_LINK_LIBRARIES(Core Weave Dirac SU3 QCD)
