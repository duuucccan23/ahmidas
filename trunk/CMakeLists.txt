CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Ahmidas CXX)

# Options for setting the type of build
OPTION(USE_MPI "Build Ahmidas for use on parallel systems, using MPI for parallelism." OFF)
OPTION(USE_OPENMP "Build Ahmidas for use on shared memory parallel systems." OFF)
OPTION(USE_SCALAR "Build Ahmidas for use on scalar systems [default]." ON)

IF(USE_MPI)
FIND_PACKAGE(MPI REQUIRED)
ENDIF(USE_MPI)

ENABLE_TESTING()

INCLUDE_DIRECTORIES(BEFORE SYSTEM src)
SET(CMAKE_CXX_FLAGS "-O3 -pg -g -Wall -pedantic")

IF(USE_SCALAR)
INCLUDE_DIRECTORIES(BEFORE SYSTEM arch/scalar)
ADD_SUBDIRECTORY(arch/scalar)
ENDIF(USE_SCALAR)

IF(USE_MPI)
FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(BEFORE SYSTEM arch/mpi)
ADD_SUBDIRECTORY(arch/mpi)
ENDIF(USE_MPI)

IF(USE_OPENMP)
INCLUDE_DIRECTORIES(BEFORE SYSTEM arch/openmp)
ADD_SUBDIRECTORY(arch/openmp)
ENDIF(USE_OPENMP)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(example)

ADD_EXECUTABLE(main src/main.cpp)
TARGET_LINK_LIBRARIES(main Weave SU3 Lime QCD ILDG MILC Base)

ADD_EXECUTABLE(simon_test src/simon_test.cpp)
TARGET_LINK_LIBRARIES(simon_test SU3 Lime QCD ILDG Base Core Weave)


