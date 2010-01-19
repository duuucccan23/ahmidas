SET(Sources_Meson Meson/Meson.cpp Meson.h)
ADD_LIBRARY(Contract ${Sources_Meson})
SET_TARGET_PROPERTIES(Contract PROPERTIES LINKER_LANGUAGE CXX)
