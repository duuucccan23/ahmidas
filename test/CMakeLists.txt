ADD_EXECUTABLE(endian endian.cpp)
ADD_TEST(endianness endian)

ADD_EXECUTABLE(reunit1 reunit1.cpp)
TARGET_LINK_LIBRARIES(reunit1 SU3)
ADD_TEST(reunitarization_1 reunit1)

ADD_EXECUTABLE(reunit2 reunit2.cpp)
TARGET_LINK_LIBRARIES(reunit2 SU3)
ADD_TEST(reunitarization_2 reunit2)

ADD_EXECUTABLE(limereader limereader.cpp)
TARGET_LINK_LIBRARIES(limereader Lime)
ADD_TEST(limereader limereader conf.2448)
