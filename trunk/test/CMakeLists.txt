ADD_EXECUTABLE(endian endian.cpp)
ADD_TEST(endianness endian)

ADD_EXECUTABLE(reunit reunit.cpp)
TARGET_LINK_LIBRARIES(reunit SU3)
ADD_TEST(reunitarization_stress_test reunit)

ADD_EXECUTABLE(weave weave.cpp)
TARGET_LINK_LIBRARIES(weave Base)
ADD_TEST(weave_initialization weave)

ADD_EXECUTABLE(limereader limereader.cpp)
TARGET_LINK_LIBRARIES(limereader Lime)
ADD_TEST(limereader_read_in limereader ../../test/conf.88)

ADD_EXECUTABLE(fieldcomponent fieldcomponent.cpp)
TARGET_LINK_LIBRARIES(fieldcomponent Core)
ADD_TEST(fieldcomponent_setup fieldcomponent)

ADD_EXECUTABLE(random random.cpp)
TARGET_LINK_LIBRARIES(random Base)
ADD_TEST(random_performance random)

IF(USE_MPI)
ADD_EXECUTABLE(MPIstartstop MPIstartstop.cpp)
ADD_TEST(MPIstarting_stopping MPIstartstop)
ENDIF(USE_MPI)
