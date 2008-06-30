CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
PROJECT(Ahmidas)

FIND_PACKAGE(MPI REQUIRED)

INCLUDE_DIRECTORIES(BEFORE SYSTEM ../src)

ADD_SUBDIRECTORY(Core)
ADD_SUBDIRECTORY(Dirac)
ADD_SUBDIRECTORY(Lime)
ADD_SUBDIRECTORY(QCD)
ADD_SUBDIRECTORY(SU3)


ADD_EXECUTABLE(main main.cpp)
TARGET_LINK_LIBRARIES(main Lime SU3)

ENABLE_TESTING()
ADD_EXECUTABLE(../../test/endian ../../test/endian.cpp)
ADD_TEST(endianness ../../test/endian)

ADD_EXECUTABLE(../../test/reunit1 ../../test/reunit1.cpp)
TARGET_LINK_LIBRARIES(../../test/reunit1 SU3)
ADD_TEST(reunitarization_1 ../../test/reunit1)

ADD_EXECUTABLE(../../test/reunit2 ../../test/reunit2.cpp)
TARGET_LINK_LIBRARIES(../../test/reunit2 SU3)
ADD_TEST(reunitarization_2 ../../test/reunit2)
