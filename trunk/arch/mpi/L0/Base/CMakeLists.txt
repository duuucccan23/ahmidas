SET(Sources_Weave Weave.h Weave/Weave_Weave_a.cpp)
ADD_LIBRARY(Weave ${Sources_Weave})
TARGET_LINK_LIBRARIES(Weave Base Grid)

SET(Sources_Grid Grid.h Grid/Grid_rank.cpp Grid/Grid_initContiguousBlock.cpp Grid/Grid_greatestCommonDivisor.cpp)
ADD_LIBRARY(Grid ${Sources_Grid})
TARGET_LINK_LIBRARIES(Grid Base)

ADD_SUBDIRECTORY(IO)