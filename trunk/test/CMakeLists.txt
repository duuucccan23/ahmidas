ADD_EXECUTABLE(endian endian.cpp)
ADD_TEST(endianness endian)

ADD_EXECUTABLE(reunit1 reunit1.cpp)
TARGET_LINK_LIBRARIES(reunit1 SU3)
ADD_TEST(reunitarization_1 reunit1)

ADD_EXECUTABLE(reunit2 reunit2.cpp)
TARGET_LINK_LIBRARIES(reunit2 SU3)
ADD_TEST(reunitarization_2 reunit2)

ADD_EXECUTABLE(weave weave.cpp)
TARGET_LINK_LIBRARIES(weave)
ADD_TEST(weave weave)

ADD_EXECUTABLE(limereader limereader.cpp)
TARGET_LINK_LIBRARIES(limereader Lime)
ADD_TEST(limereader limereader ../../test/conf.88)

IF(USE_MPI) # THIS SHOULD NOT BE NECESSARY!
ADD_EXECUTABLE(fieldcomponent fieldcomponent.cpp)
TARGET_LINK_LIBRARIES(fieldcomponent Core)
ADD_TEST(fieldcomponent fieldcomponent)
ENDIF(USE_MPI)

IF(USE_MPI)
ADD_EXECUTABLE(MPIstartstop MPIstartstop.cpp)
ADD_TEST(MPIstarting_stopping MPIstartstop)
ENDIF(USE_MPI)
