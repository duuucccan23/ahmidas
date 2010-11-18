ADD_SUBDIRECTORY(Smear)

SET(Sources_Tool Tool/Tool_localTrace.cpp Tool/Tool_reunitarize.cpp Tool/Tool_fixCoulombGauge.cpp
                 Tool/Tool_temporalPlaquette.cpp Tool/Tool_localRealTrace.cpp Tool/Tool_realtr.cpp
                 Tool/Tool_spatialPlaquette.cpp Tool/Tool_tr.cpp Tool/ScidacChecksum/ScidacChecksum_createTable.cpp
                 Tool/Tool_printLightMesonCorrelator.cpp Tool/Tool_spatialDownPlaquette.cpp 
                 Tool/Tool_temporalDownPlaquette.cpp
                 Tool/Tool_convertIntToOperator.cpp)
ADD_LIBRARY(Tool ${Sources_Tool})
TARGET_LINK_LIBRARIES(Tool SU3 Base Core Path Weave Print)

SET(Sources_Path Path/Path_square.cpp Path/Path_staple.cpp Path/Path_step.cpp)
ADD_LIBRARY(Path ${Sources_Path})
TARGET_LINK_LIBRARIES(Path QCD SU3 Base Weave)
