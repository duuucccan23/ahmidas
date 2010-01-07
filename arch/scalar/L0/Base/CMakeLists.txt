SET(Sources_Weave Weave/Weave_Weave_a.cpp Weave/Weave_Weave_b.cpp Weave/Weave_operatorassign.cpp)
ADD_LIBRARY(Weave ${Sources_Weave})
TARGET_LINK_LIBRARIES(Weave Base)
