ADD_SUBDIRECTORY(IO)

SET(Sources_ILDG   IO/ILDGinfo_ILDGinfo.cpp)
SET(Sources_Scidac IO/ScidacInfo_ScidacInfo.cpp)
SET(Sources_IO     IO/load.cpp IO/save.cpp)

IF (${LEMON_FOUND})
  SET(Sources_IO_Lemon  IO/loadPropagator_a.cpp IO/loadPropagator_b.cpp IO/loadScidacInfo.cpp IO/savePropagator.cpp IO/savePropagatorInfo.cpp IO/savePropagatorType.cpp IO/ScidacInfo_parseScidacInfo.cpp)
ELSE (${LEMON_FOUND})
  SET(Sources_IO_Lemon "")
ENDIF (${LEMON_FOUND})

ADD_LIBRARY(IO ${Sources_IO_Lemon} ${Sources_ILDG} ${Sources_Scidac} ${Sources_IO})
TARGET_LINK_LIBRARIES(IO Tool Lime)

