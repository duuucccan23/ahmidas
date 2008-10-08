SET(Sources_Base Base.h)

IF(USE_MPI)
SET(Sources_Grid Grid.h)
SET(Sources_Unifold Unifold.h)
ENDIF(USE_MPI)

ADD_LIBRARY(Base ${Sources_Base} ${Sources_Grid} ${Sources_Unifold})
SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)

ADD_SUBDIRECTORY(IO)