ADD_EXECUTABLE(generate_sequential_sources_nucleon generate_sequential_sources_nucleon.cpp)
TARGET_LINK_LIBRARIES(generate_sequential_sources_nucleon SU3 QCD IO Base Core Dirac Weave Contract Input Smear Ahmidas)

ADD_EXECUTABLE(contract_nucleon3point contract_nucleon3point.cpp)
TARGET_LINK_LIBRARIES(contract_nucleon3point SU3 QCD IO Base Core Dirac Weave Contract Input Smear Ahmidas lemon)

ADD_EXECUTABLE(contract_nucleon3point_stochastic contract_nucleon3point_stochastic.cpp)
IF(USE_MPI)
TARGET_LINK_LIBRARIES(contract_nucleon3point_stochastic SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas lemon)
ELSE(USE_MPI)
TARGET_LINK_LIBRARIES(contract_nucleon3point_stochastic SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas)
ENDIF(USE_MPI)

ADD_EXECUTABLE(contract_baryon2point contract_baryon2point.cpp)
IF(USE_MPI)
TARGET_LINK_LIBRARIES(contract_baryon2point SU3 QCD IO Base Core Dirac Weave Contract Input Smear Ahmidas lemon)
ELSE(USE_MPI)
TARGET_LINK_LIBRARIES(contract_baryon2point SU3 QCD IO Base Core Dirac Weave Contract Input Smear Ahmidas)
ENDIF(USE_MPI)

ADD_EXECUTABLE(contract_meson_2pts contract_meson_2pts.cpp)
TARGET_LINK_LIBRARIES(contract_meson_2pts SU3 QCD IO Base Core Dirac Contract Input Ahmidas Print)

ADD_EXECUTABLE(contract_meson_stochastic_2pts contract_meson_stochastic_2pts.cpp)
TARGET_LINK_LIBRARIES(contract_meson_stochastic_2pts SU3 QCD IO Base Core Dirac Contract Input Ahmidas Print)

ADD_EXECUTABLE(contract_meson_ultrastochastic_2pts contract_meson_ultrastochastic_2pts.cpp)
TARGET_LINK_LIBRARIES(contract_meson_ultrastochastic_2pts SU3 QCD IO Base Core Dirac Contract Input Ahmidas Print)


ADD_EXECUTABLE(contract_two_lines contract_two_lines.cpp)
TARGET_LINK_LIBRARIES(contract_two_lines SU3 QCD IO Base Core Dirac Contract Input Ahmidas Print)

ADD_EXECUTABLE(contract_two_stochastic_lines contract_two_stochastic_lines.cpp)
TARGET_LINK_LIBRARIES(contract_two_stochastic_lines SU3 QCD IO Base Core Dirac Contract Input Ahmidas Print)

ADD_EXECUTABLE(contract_two_ultrastochastic_lines contract_two_ultrastochastic_lines.cpp)
TARGET_LINK_LIBRARIES(contract_two_ultrastochastic_lines SU3 QCD IO Base Core Dirac Contract Input Ahmidas Print)


ADD_EXECUTABLE(calculate_BK calculate_BK.cpp)
TARGET_LINK_LIBRARIES(calculate_BK SU3 QCD IO Base Core Dirac Weave Contract Input Ahmidas Print)


ADD_EXECUTABLE(generate_meson_sequential_source generate_meson_sequential_source.cpp)
TARGET_LINK_LIBRARIES(generate_meson_sequential_source IO Input Core Ahmidas)

ADD_EXECUTABLE(generate_meson_sequential_stochastic_source generate_meson_sequential_stochastic_source.cpp)
TARGET_LINK_LIBRARIES(generate_meson_sequential_stochastic_source IO Input Core Ahmidas)

ADD_EXECUTABLE(generate_meson_sequential_ultrastochastic_source generate_meson_sequential_ultrastochastic_source.cpp)
TARGET_LINK_LIBRARIES(generate_meson_sequential_ultrastochastic_source IO Input Core Ahmidas)


ADD_EXECUTABLE(finalize_meson_sequential_propagator finalize_meson_sequential_propagator.cpp)
TARGET_LINK_LIBRARIES(finalize_meson_sequential_propagator IO Input Core Ahmidas)

ADD_EXECUTABLE(finalize_meson_sequential_stochastic_propagator finalize_meson_sequential_stochastic_propagator.cpp)
TARGET_LINK_LIBRARIES(finalize_meson_sequential_stochastic_propagator IO Input Core Ahmidas)

ADD_EXECUTABLE(finalize_meson_sequential_ultrastochastic_propagator finalize_meson_sequential_ultrastochastic_propagator.cpp)
TARGET_LINK_LIBRARIES(finalize_meson_sequential_ultrastochastic_propagator IO Input Core Ahmidas)

ADD_EXECUTABLE(contract_disconnected contract_disconnected.cpp)
TARGET_LINK_LIBRARIES(contract_disconnected Print Ahmidas SU3 QCD IO Base Core Dirac
			Weave Contract Input Smear)

ADD_EXECUTABLE(contract_disconnected_vvSum contract_disconnected_vvSum.cpp)
TARGET_LINK_LIBRARIES(contract_disconnected_vvSum Print Ahmidas SU3 QCD IO Base Core Dirac
			Weave Contract Input Smear)

ADD_EXECUTABLE(contract_disconnected_heavylight contract_disconnected_heavylight.cpp)
TARGET_LINK_LIBRARIES(contract_disconnected_heavylight Print Ahmidas SU3 QCD IO Base Core Dirac
			Weave Contract Input Smear)


