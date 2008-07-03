ADD_EXECUTABLE(endian endian.cpp)
ADD_TEST(endianness endian)

ADD_EXECUTABLE(grid grid.cpp)
TARGET_LINK_LIBRARIES(grid Core)
ADD_TEST(grid grid)

ADD_EXECUTABLE(reunit1 reunit1.cpp)
TARGET_LINK_LIBRARIES(reunit1 SU3)
ADD_TEST(reunitarization_1 reunit1)

ADD_EXECUTABLE(reunit2 reunit2.cpp)
TARGET_LINK_LIBRARIES(reunit2 SU3)
ADD_TEST(reunitarization_2 reunit2)

ADD_EXECUTABLE(limereader limereader.cpp)
TARGET_LINK_LIBRARIES(limereader Lime)
ADD_TEST(limereader limereader ../../test/conf.2448 ildg-binary-data)

ADD_EXECUTABLE(fieldcomponent fieldcomponent.cpp)
TARGET_LINK_LIBRARIES(fieldcomponent Core)
ADD_TEST(fieldcomponent fieldcomponent)

ADD_EXECUTABLE(fieldIncreaseIndex fieldIncreaseIndex.cpp)
TARGET_LINK_LIBRARIES(fieldIncreaseIndex Core)
ADD_TEST(fieldIncreaseIndex fieldIncreaseIndex)
