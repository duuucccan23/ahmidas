SET(Sources_Component Component.h)
SET(Sources_Correlator Correlator.h)
SET(Sources_Field Field.h)
SET(Sources_Propagator Propagator.h Propagator/StochasticSource.cpp)
ADD_LIBRARY(Core ${Sources_Component} ${Sources_Correlator} ${Sources_Field} ${Sources_Propagator})
SET_TARGET_PROPERTIES(Core PROPERTIES LINKER_LANGUAGE CXX)

