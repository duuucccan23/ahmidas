SET(Sources_Weave Weave/Weave.cpp Weave.h Weave/Weave_instance.cpp Weave/Weave_Weave.cpp)
ADD_LIBRARY(Base ${Sources_Weave})
SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)