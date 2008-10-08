CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Ahmidas CXX)

# Options for setting the type of build
OPTION(USE_MPI "Build Ahmidas for use on parallel systems, using MPI for parallelism." OFF)
OPTION(USE_OPENMP "Build Ahmidas for use on shared memory parallel systems." OFF)
OPTION(USE_SCALAR "Build Ahmidas for use on scalar systems [default]." ON)

IF (USE_MPI)
FIND_PACKAGE(MPI REQUIRED)
ENDIF (USE_MPI)

ENABLE_TESTING()

INCLUDE_DIRECTORIES(BEFORE SYSTEM src)

IF (USE_SCALAR)
INCLUDE_DIRECTORIES(BEFORE SYSTEM arch)
ENDIF (USE_SCALAR)

IF (USE_MPI)
FIND_PACKAGE(MPI REQUIRED)
ENDIF (USE_MPI)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)

ADD_EXECUTABLE(main src/main.cpp)
TARGET_LINK_LIBRARIES(main)
