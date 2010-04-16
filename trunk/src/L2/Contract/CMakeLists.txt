SET(Sources_Meson Meson.h Meson/Meson.cpp)
SET(Sources_Baryon Baryon.h
                   Baryon/Proton.cpp Baryon/Baryon_proton_threepoint_naive.cpp
                   Baryon/Baryon_create_sequential_sources.cpp)
ADD_LIBRARY(Contract ${Sources_Meson} ${Sources_Baryon})
SET_TARGET_PROPERTIES(Contract PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(Contract Core Dirac QCD)