SET(Sources_Meson Meson.h Meson/Meson.cpp)
SET(Sources_Baryon Baryon.h
                   Baryon/Baryon_proton_twopoint.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_a.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_b.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_c.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_d.cpp
                   Baryon/Baryon_create_sequential_source_proton_fixed_insertion_timeslice.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_a.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_b.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_c.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_d.cpp
                   Baryon/Baryon_proton_threepoint_sequential.cpp
                   Baryon/Baryon_proton_threepoint_stochastic_naive.cpp
                   Baryon/Baryon_proton_threepoint_stochastic.cpp)
ADD_LIBRARY(Contract ${Sources_Meson} ${Sources_Baryon})
SET_TARGET_PROPERTIES(Contract PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(Contract Core Dirac QCD)