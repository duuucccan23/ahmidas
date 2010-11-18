ADD_EXECUTABLE(endian endian.cpp)
TARGET_LINK_LIBRARIES(endian Base Ahmidas Print)
ADD_TEST(endianness endian)

ADD_EXECUTABLE(isnan isnan.cpp)
TARGET_LINK_LIBRARIES(isnan Base Ahmidas Print)
ADD_TEST(isnan isnan)

ADD_EXECUTABLE(reunit reunit.cpp)
TARGET_LINK_LIBRARIES(reunit SU3 Ahmidas Print Debug)
ADD_TEST(reunitarization_stress_test reunit)

ADD_EXECUTABLE(weave weave.cpp)
TARGET_LINK_LIBRARIES(weave Base Ahmidas Weave Debug)
ADD_TEST(weave_initialization weave)

ADD_EXECUTABLE(field_shift field_shift.cpp)
TARGET_LINK_LIBRARIES(field_shift SU3 QCD Ahmidas Base Core Weave)
ADD_TEST(field_shift_test field_shift)

ADD_EXECUTABLE(limereader limereader.cpp)
TARGET_LINK_LIBRARIES(limereader Lime Base Weave)
ADD_TEST(limereader_read_in limereader ../../test/conf.88)

ADD_EXECUTABLE(readFile readFile.cpp)
TARGET_LINK_LIBRARIES(readFile Input Ahmidas)
ADD_TEST(read_input_file readFile)

ADD_EXECUTABLE(loadMILC loadMILC.cpp)
TARGET_LINK_LIBRARIES(loadMILC Base IO SU3 Ahmidas)
ADD_TEST(load_MILC_lattice loadMILC)

ADD_EXECUTABLE(random random.cpp)
TARGET_LINK_LIBRARIES(random Base Print Ahmidas)
ADD_TEST(random_performance random)

ADD_EXECUTABLE(config config.cpp)
TARGET_LINK_LIBRARIES(config Base IO SU3 Print Ahmidas)
ADD_TEST(basic_config_stuff config)

ADD_EXECUTABLE(light_mesons light_mesons.cpp)
TARGET_LINK_LIBRARIES(light_mesons SU3 QCD IO Base Core Dirac Weave Contract Ahmidas)
ADD_TEST(light_mesons_naive light_mesons)

ADD_EXECUTABLE(light_mesons_one_end_trick light_mesons_one_end_trick.cpp)
TARGET_LINK_LIBRARIES(light_mesons_one_end_trick SU3 QCD IO Base Core Dirac Weave Contract Ahmidas)
ADD_TEST(light_mesons_one_end light_mesons_one_end_trick)

ADD_EXECUTABLE(proton proton.cpp)
TARGET_LINK_LIBRARIES(proton SU3 QCD IO Base Core Weave Contract Ahmidas)
ADD_TEST(proton_two_point proton)

ADD_EXECUTABLE(scidacchecksum scidacchecksum.cpp)
TARGET_LINK_LIBRARIES(scidacchecksum QCD IO Print Ahmidas)
ADD_TEST(scidac_checksum scidacchecksum)

ADD_EXECUTABLE(sequential_sources sequential_sources.cpp)
TARGET_LINK_LIBRARIES(sequential_sources SU3 QCD IO Base Core Weave Contract Ahmidas)
ADD_TEST(generate_sequential_source sequential_sources)

ADD_EXECUTABLE(smearing smearing.cpp)
TARGET_LINK_LIBRARIES(smearing SU3 QCD IO Base Core Weave Smear Ahmidas)
ADD_TEST(smearing_test smearing)

ADD_EXECUTABLE(dirac_operator dirac_operator.cpp)
TARGET_LINK_LIBRARIES(dirac_operator SU3 QCD IO Base Core Weave Ahmidas)
ADD_TEST(dirac_operator_test dirac_operator)

IF(EIGEN3_FOUND)
  ADD_EXECUTABLE(reunit_algo_comp reunit_algo_comp.cpp)
  TARGET_LINK_LIBRARIES(reunit_algo_comp SU3 Ahmidas)
  ADD_TEST(reunitarization_algorithm_comparison reunit_algo_comp)
ENDIF(EIGEN3_FOUND)

# IO test for parallel IO mainly
# this is not made a test yet because it is unclear
# if the lattice can be distributed to the number of processes specified
ADD_EXECUTABLE(testIO testIO.cpp)
TARGET_LINK_LIBRARIES(testIO IO Base Core Weave Ahmidas)

# ADD_EXECUTABLE(proton_3point_sequential proton_3point_sequential.cpp)
# TARGET_LINK_LIBRARIES(proton_3point_sequential SU3 QCD IO Base Core Weave Contract)
# ADD_TEST(proton_threepoint_sequential 

# IF(USE_MPI)
# ADD_EXECUTABLE(MPIstartstop MPIstartstop.cpp)
# ADD_TEST(MPIstarting_stopping MPIstartstop Ahmidas)
# ENDIF(USE_MPI)
