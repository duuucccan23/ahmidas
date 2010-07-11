ADD_EXECUTABLE(confrw confrw.cpp)
TARGET_LINK_LIBRARIES(confrw Base IO SU3 Tool)

ADD_EXECUTABLE(smearAPE smearAPE.cpp)
TARGET_LINK_LIBRARIES(smearAPE Tool Path Base IO SU3 Weave Smear)

ADD_EXECUTABLE(pion_contr pion_contraction.cpp)
TARGET_LINK_LIBRARIES(pion_contr Base Core IO Contract Input)

ADD_EXECUTABLE(gamma gamma.cpp)
TARGET_LINK_LIBRARIES(gamma Dirac)

ADD_EXECUTABLE(generate_sequential_sources generate_sequential_sources.cpp)
TARGET_LINK_LIBRARIES(generate_sequential_sources SU3 IO Base Core Dirac Weave Contract Input Smear)

ADD_EXECUTABLE(generate_sequential_sources_alternative generate_sequential_sources_alternative.cpp)
TARGET_LINK_LIBRARIES(generate_sequential_sources_alternative SU3 IO Base Core Dirac Weave Contract Input Smear)

ADD_EXECUTABLE(contract_3point_sequential contract_3point_sequential.cpp)
TARGET_LINK_LIBRARIES(contract_3point_sequential SU3 IO Base Core Dirac Weave Contract Input Smear)

ADD_EXECUTABLE(contract_3point_sequential_alternative contract_3point_sequential_alternative.cpp)
TARGET_LINK_LIBRARIES(contract_3point_sequential_alternative SU3 IO Base Core Dirac Weave Contract Input Smear)

ADD_EXECUTABLE(prova prova.cpp)
TARGET_LINK_LIBRARIES(prova Base Core IO Input)

ADD_EXECUTABLE(prova_contraction prova_contraction.cpp)
TARGET_LINK_LIBRARIES(prova_contraction Base Core IO Input)

ADD_EXECUTABLE(prova_sequential_contraction prova_sequential_contraction.cpp)
TARGET_LINK_LIBRARIES(prova_sequential_contraction Base Core IO Input)

ADD_EXECUTABLE(generate_stochastic_source generate_stochastic_source.cpp)
TARGET_LINK_LIBRARIES(generate_stochastic_source SU3 IO Base Core Dirac Weave Input Smear)

ADD_EXECUTABLE(generate_point_source generate_point_source.cpp)
TARGET_LINK_LIBRARIES(generate_point_source SU3 IO Base Core Dirac Weave Input Smear)