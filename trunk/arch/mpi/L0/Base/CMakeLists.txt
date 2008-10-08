SET(Sources_Base Base.h)
SET(Sources_Grid Grid.h)
SET(Sources_Unifold Unifold.h)

ADD_LIBRARY(Base ${Sources_Base} ${Sources_Grid} ${Sources_Unifold})
SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)

ADD_SUBDIRECTORY(IO)