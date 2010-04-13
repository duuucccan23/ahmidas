SET(Sources_Correlator Correlator/Correlator_sumOverSpatialVolume.cpp Correlator/Correlator_projection.cpp)
SET(Sources_Propagator Propagator/StochasticSource.cpp
    Propagator/Propagator_changeBoundaryConditions.cpp
    Propagator/Propagator_gaugeTransformation.cpp
    Propagator/Propagator_construct_baryon_with_operator_insertion.cpp)
ADD_LIBRARY(Core ${Sources_Correlator} ${Sources_Propagator})
TARGET_LINK_LIBRARIES(Core Weave Dirac SU3 QCD)
