SET(Sources_Component Component.h)
SET(Sources_Correlator Correlator.h Correlator/Correlator_sumOverSpatialVolume.cpp Correlator/Correlator_projection.cpp)
SET(Sources_Field Field.h)
SET(Sources_Propagator Propagator.h Propagator/StochasticSource.cpp
    Propagator/Propagator_changeBoundaryConditions.cpp
    Propagator/Propagator_gaugeTransformation.cpp
    Propagator/Propagator_construct_baryon_with_operator_insertion.cpp)
ADD_LIBRARY(Core ${Sources_Component} ${Sources_Correlator} ${Sources_Field} ${Sources_Propagator})
TARGET_LINK_LIBRARIES(Core Weave Dirac SU3 QCD)
