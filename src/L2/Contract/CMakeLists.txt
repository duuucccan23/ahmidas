SET(Sources_Meson Meson.h
                  Meson/Meson_light_meson_twopoint_alternative.cpp
                  Meson/Meson_light_meson_twopoint_stochastic_a.cpp
                  Meson/Meson_light_meson_twopoint_stochastic_b.cpp
                  Meson/Meson_light_meson_twopoint_ultrastochastic.cpp
		  Meson/Meson_BK_threepoint_disconnected_stochastic.cpp
		  Meson/Meson_BK_threepoint_connected_stochastic.cpp)
SET(Sources_Baryon Baryon.h
                   Baryon/Baryon_proton_twopoint.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_a.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_b.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_c.cpp
                   Baryon/Baryon_create_sequential_source_proton_d_d.cpp
                   Baryon/Baryon_create_sequential_source_proton_fixed_insertion_timeslice.cpp
                   Baryon/Baryon_create_sequential_source_proton_fixed_insertion_timeslice_derivative_operator.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_a.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_b.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_c.cpp
                   Baryon/Baryon_create_sequential_source_proton_u_d.cpp
                   Baryon/Baryon_proton_threepoint_sequential.cpp
                   Baryon/Baryon_proton_threepoint_stochastic_naive.cpp
                   Baryon/Baryon_proton_threepoint_stochastic.cpp)
SET(Sources_Disconnected Disconnected.h
                  Disconnected/Disconnected_compute_loop_twist2_operator_a.cpp
                  Disconnected/Disconnected_compute_loop_twist2_operator_b.cpp
				  Disconnected/Disconnected_compute_loop_conserved_current.cpp
				  Disconnected/Disconnected_compute_loop_a.cpp
                  Disconnected/Disconnected_compute_loop_b.cpp
                  Disconnected/Disconnected_compute_loop_c.cpp)

ADD_LIBRARY(Contract ${Sources_Meson} ${Sources_Baryon} ${Sources_Disconnected})
SET_TARGET_PROPERTIES(Contract PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(Contract Core Dirac QCD)
