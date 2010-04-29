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

ADD_EXECUTABLE(readFile readFile.cpp)
TARGET_LINK_LIBRARIES(readFile Input)
ADD_TEST(read_input_file readFile)

ADD_EXECUTABLE(fieldcomponent fieldcomponent.cpp)
TARGET_LINK_LIBRARIES(fieldcomponent Weave)
ADD_TEST(fieldcomponent_setup fieldcomponent)

ADD_EXECUTABLE(loadMILC loadMILC.cpp)
TARGET_LINK_LIBRARIES(loadMILC Base IO SU3)
ADD_TEST(load_MILC_lattice loadMILC)

ADD_EXECUTABLE(random random.cpp)
TARGET_LINK_LIBRARIES(random Base)
ADD_TEST(random_performance random)

ADD_EXECUTABLE(config config.cpp)
TARGET_LINK_LIBRARIES(config Base IO SU3)
ADD_TEST(basic_config_stuff config)

ADD_EXECUTABLE(light_mesons light_mesons.cpp)
TARGET_LINK_LIBRARIES(light_mesons SU3 QCD IO Base Core Dirac Weave Contract)
ADD_TEST(light_mesons_naive light_mesons)

ADD_EXECUTABLE(light_mesons_one_end_trick light_mesons_one_end_trick.cpp)
TARGET_LINK_LIBRARIES(light_mesons_one_end_trick SU3 QCD IO Base Core Dirac Weave Contract)
ADD_TEST(light_mesons_one_end light_mesons_one_end_trick)

ADD_EXECUTABLE(proton proton.cpp)
TARGET_LINK_LIBRARIES(proton SU3 QCD IO Base Core Weave Contract)
ADD_TEST(proton_two_point proton)

ADD_EXECUTABLE(sequential_sources sequential_sources.cpp)
TARGET_LINK_LIBRARIES(sequential_sources SU3 QCD IO Base Core Weave Contract)
ADD_TEST(generate_sequential_source sequential_sources)

# IF(USE_MPI)
# ADD_EXECUTABLE(MPIstartstop MPIstartstop.cpp)
# ADD_TEST(MPIstarting_stopping MPIstartstop)
# ENDIF(USE_MPI)
