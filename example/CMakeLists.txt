ADD_EXECUTABLE(confrw confrw.cpp)
TARGET_LINK_LIBRARIES(confrw Base ILDG Lime SU3 Tool)

ADD_EXECUTABLE(smearAPE smearAPE.cpp)
TARGET_LINK_LIBRARIES(smearAPE Tool Path Base ILDG Lime SU3 Weave Smear)
