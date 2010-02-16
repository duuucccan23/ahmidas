SET(Sources_Meson Meson.h Meson/Meson.cpp)
ADD_LIBRARY(Contract ${Sources_Meson})
SET_TARGET_PROPERTIES(Contract PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(Contract Core Dirac)