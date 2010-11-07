SET(Sources_Correlator Correlator/Correlator.cpp
                       Correlator/Correlator_projection.cpp
                       Correlator/Correlator_cout_operator_lshift.cpp)
SET(Sources_Propagator Propagator/StochasticPropagator.cpp
                       Propagator/StochasticSource.cpp
                       Propagator/StochasticPropagator_rotateToPhysicalBasis.cpp
                       Propagator/Propagator_applyDiracOperator_a.cpp
                       Propagator/Propagator_applyDiracOperator_b.cpp
                       Propagator/Propagator_applyDiracOperator_c.cpp
                       Propagator/Propagator_changeBoundaryConditions.cpp
                       Propagator/Propagator_construct_baryon.cpp
                       Propagator/Propagator_construct_baryon_with_operator_insertion.cpp
                       Propagator/Propagator_contract_a.cpp
                       Propagator/Propagator_destroy.cpp
                       Propagator/Propagator_diff.cpp
                       Propagator/Propagator_gaugeTransformation.cpp
                       Propagator/Propagator_isolate.cpp
                       Propagator/Propagator_multiplyOperator.cpp
                       Propagator/Propagator_norm.cpp
                       Propagator/Propagator_normq.cpp
                       Propagator/Propagator_reconstruct_doublet.cpp
                       Propagator/Propagator_revert.cpp
                       Propagator/Propagator_rotateToPhysicalBasis.cpp
                       Propagator/Propagator_select_timeslice.cpp
                       Propagator/Propagator_shift.cpp
                       Propagator/Propagator_smearJacobi.cpp
                       Propagator/Propagator_trace.cpp
                       Propagator/Propagator.operators.cpp)
ADD_LIBRARY(Core ${Sources_Correlator} ${Sources_Propagator})
TARGET_LINK_LIBRARIES(Core Weave Dirac SU3 QCD)
