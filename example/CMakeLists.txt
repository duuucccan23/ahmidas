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
