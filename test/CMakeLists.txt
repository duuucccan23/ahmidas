ADD_EXECUTABLE(endian endian.cpp)
TARGET_LINK_LIBRARIES(Core)
ADD_TEST(endianness endian)

IF(USE_MPI) # THIS SHOULD NOT BE NECESSARY!
ADD_EXECUTABLE(grid grid.cpp)
TARGET_LINK_LIBRARIES(grid Core)
ADD_TEST(grid grid)
ENDIF(USE_MPI)

ADD_EXECUTABLE(reunit1 reunit1.cpp)
TARGET_LINK_LIBRARIES(reunit1 SU3)
ADD_TEST(reunitarization_1 reunit1)

ADD_EXECUTABLE(reunit2 reunit2.cpp)
TARGET_LINK_LIBRARIES(reunit2 SU3)
ADD_TEST(reunitarization_2 reunit2)

IF(USE_MPI) # THIS SHOULD NOT BE NECESSARY!
ADD_EXECUTABLE(limereader limereader.cpp)
TARGET_LINK_LIBRARIES(limereader Lime)
ADD_TEST(limereader limereader ../../test/conf.88)
ENDIF(USE_MPI)

IF(USE_MPI) # THIS SHOULD NOT BE NECESSARY!
ADD_EXECUTABLE(fieldcomponent fieldcomponent.cpp)
TARGET_LINK_LIBRARIES(fieldcomponent Core)
ADD_TEST(fieldcomponent fieldcomponent)
ENDIF(USE_MPI)

IF(USE_MPI) # THIS SHOULD NOT BE NECESSARY!
ADD_EXECUTABLE(fieldIncreaseIndex fieldIncreaseIndex.cpp)
TARGET_LINK_LIBRARIES(fieldIncreaseIndex Core)
ADD_TEST(fieldincreaseindex fieldIncreaseIndex)
ENDIF(USE_MPI)

IF(USE_MPI)
ADD_EXECUTABLE(MPIstartstop MPIstartstop.cpp)
ADD_TEST(MPIstarting_stopping MPIstartstop)
ENDIF(USE_MPI)
