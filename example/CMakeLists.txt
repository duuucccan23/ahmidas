ADD_EXECUTABLE(confrw confrw.cpp)
TARGET_LINK_LIBRARIES(confrw Base IO SU3 Tool Ahmidas)

ADD_EXECUTABLE(smearAPE smearAPE.cpp)
TARGET_LINK_LIBRARIES(smearAPE Tool Path Base IO SU3 Weave Smear Ahmidas)

ADD_EXECUTABLE(pion_contr pion_contraction.cpp)
TARGET_LINK_LIBRARIES(pion_contr Base Core IO Contract Input Ahmidas)

ADD_EXECUTABLE(gamma gamma.cpp)
TARGET_LINK_LIBRARIES(gamma Dirac Ahmidas)

ADD_EXECUTABLE(generate_sequential_sources generate_sequential_sources.cpp)
TARGET_LINK_LIBRARIES(generate_sequential_sources SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas)

ADD_EXECUTABLE(generate_sequential_sources_alternative generate_sequential_sources_alternative.cpp)
TARGET_LINK_LIBRARIES(generate_sequential_sources_alternative SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas)

ADD_EXECUTABLE(contract_3point_sequential contract_3point_sequential.cpp)
TARGET_LINK_LIBRARIES(contract_3point_sequential SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas lemon)

ADD_EXECUTABLE(contract_2point_proton contract_2point_proton.cpp)
TARGET_LINK_LIBRARIES(contract_2point_proton SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas)

ADD_EXECUTABLE(contract_3point_sequential_alternative contract_3point_sequential_alternative.cpp)
TARGET_LINK_LIBRARIES(contract_3point_sequential_alternative SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas)

ADD_EXECUTABLE(prova prova.cpp)
TARGET_LINK_LIBRARIES(prova Base Core IO Input Ahmidas)

ADD_EXECUTABLE(prova_contraction prova_contraction.cpp)
TARGET_LINK_LIBRARIES(prova_contraction Base Core IO Input Ahmidas)

ADD_EXECUTABLE(prova_sequential_contraction prova_sequential_contraction.cpp)
TARGET_LINK_LIBRARIES(prova_sequential_contraction Base Core IO Input Ahmidas) 

ADD_EXECUTABLE(generate_volume_source generate_volume_source.cpp)
TARGET_LINK_LIBRARIES(generate_volume_source SU3 IO Base Core Dirac Weave Input Smear Ahmidas)

ADD_EXECUTABLE(generate_stochastic_source generate_stochastic_source.cpp)
IF(USE_MPI)
  TARGET_LINK_LIBRARIES(generate_stochastic_source SU3 IO Base Core Dirac Weave Input Smear Ahmidas lemon)
ELSE(USE_MPI)
  TARGET_LINK_LIBRARIES(generate_stochastic_source SU3 IO Base Core Dirac Weave Input Smear Ahmidas)
ENDIF(USE_MPI)


ADD_EXECUTABLE(generate_point_source generate_point_source.cpp)
TARGET_LINK_LIBRARIES(generate_point_source SU3 IO Base Core Dirac Weave Input Smear Ahmidas lemon)

ADD_EXECUTABLE(output output.cpp)
TARGET_LINK_LIBRARIES(output Ahmidas Debug Print)

ADD_EXECUTABLE(convert_ildg_to_eo convert_ildg_to_eo.cpp)
TARGET_LINK_LIBRARIES(convert_ildg_to_eo IO Input Core Ahmidas)


ADD_EXECUTABLE(reconstruct_doublet reconstruct_doublet.cpp)
TARGET_LINK_LIBRARIES(reconstruct_doublet IO Input Core Ahmidas)

ADD_EXECUTABLE(reconstruct_stochastic_doublet reconstruct_stochastic_doublet.cpp)
TARGET_LINK_LIBRARIES(reconstruct_stochastic_doublet IO Input Core Ahmidas)

ADD_EXECUTABLE(reconstruct_ultrastochastic_doublet reconstruct_ultrastochastic_doublet.cpp)
TARGET_LINK_LIBRARIES(reconstruct_ultrastochastic_doublet IO Input Core Ahmidas)


ADD_EXECUTABLE(check_meson_3pts check_meson_3pts.cpp)
TARGET_LINK_LIBRARIES(check_meson_3pts IO Input Core Ahmidas)

ADD_EXECUTABLE(check_meson_stochastic_3pts check_meson_stochastic_3pts.cpp)
TARGET_LINK_LIBRARIES(check_meson_stochastic_3pts IO Input Core Ahmidas)

ADD_EXECUTABLE(check_meson_ultrastochastic_3pts check_meson_ultrastochastic_3pts.cpp)
TARGET_LINK_LIBRARIES(check_meson_ultrastochastic_3pts IO Input Core Ahmidas)


ADD_EXECUTABLE(propagator_dump propagator_dump.cpp)
TARGET_LINK_LIBRARIES(propagator_dump IO Input Core Ahmidas)
# ADD_EXECUTABLE(makeetmcprop makeetmcprop.cpp)
# TARGET_LINK_LIBRARIES(makeetmcprop IO Core Ahmidas)
