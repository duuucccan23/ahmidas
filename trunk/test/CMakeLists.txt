ADD_EXECUTABLE(endian endian.cpp)
TARGET_LINK_LIBRARIES(endian Base)
ADD_TEST(endianness endian)

ADD_EXECUTABLE(reunit reunit.cpp)
TARGET_LINK_LIBRARIES(reunit SU3)
ADD_TEST(reunitarization_stress_test reunit)

ADD_EXECUTABLE(weave weave.cpp)
TARGET_LINK_LIBRARIES(weave Base Weave)
ADD_TEST(weave_initialization weave)

ADD_EXECUTABLE(limereader limereader.cpp)
TARGET_LINK_LIBRARIES(limereader Lime Base)
ADD_TEST(limereader_read_in limereader ../../test/conf.88)

ADD_EXECUTABLE(fieldcomponent fieldcomponent.cpp)
TARGET_LINK_LIBRARIES(fieldcomponent Weave)
ADD_TEST(fieldcomponent_setup fieldcomponent)

ADD_EXECUTABLE(loadMILC loadMILC.cpp)
TARGET_LINK_LIBRARIES(loadMILC Base Tool MILC SU3)
ADD_TEST(load_MILC_lattice loadMILC)

ADD_EXECUTABLE(random random.cpp)
TARGET_LINK_LIBRARIES(random Base)
ADD_TEST(random_performance random)

ADD_EXECUTABLE(config config.cpp)
TARGET_LINK_LIBRARIES(config Base ILDG Lime SU3 Tool)
ADD_TEST(basic_config_stuff config)

#ADD_EXECUTABLE(light_mesons light_mesons.cpp)
#TARGET_LINK_LIBRARIES(light_mesons SU3 Lime QCD ILDG Scidac Base Core Dirac Weave Contract)

IF(USE_MPI)
ADD_EXECUTABLE(MPIstartstop MPIstartstop.cpp)
ADD_TEST(MPIstarting_stopping MPIstartstop)
ENDIF(USE_MPI)
