CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Ahmidas CXX)

# SET(CMAKE_CXX_COMPILER "icpc")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# Option for building the tests
OPTION(BUILD_TESTS "Build all tests packaged with Ahmidas" ON)

INCLUDE_DIRECTORIES(BEFORE SYSTEM src)
SET(CMAKE_CXX_FLAGS "-O3")

FIND_PACKAGE(Eigen3)
IF(EIGEN3_FOUND)
  MESSAGE(STATUS "Found Eigen3 header files, activating.")
  INCLUDE_DIRECTORIES(BEFORE SYSTEM ${EIGEN3_INCLUDE_DIR})  
ELSE(EIGEN3_FOUND)
  MESSAGE(STATUS "Did not find optional library Eigen3, tests will be deactivated.")
ENDIF(EIGEN3_FOUND)

# Options for setting the type of build
OPTION(USE_MPI "Build Ahmidas for use on parallel systems, using MPI for parallelism." OFF)

IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  FIND_PACKAGE(LEMON)
  IF(LEMON_FOUND)
    MESSAGE(STATUS "Lemon available, so I/O will be implemented over MPI.")
  ELSE(LEMON_FOUND)
    MESSAGE(STATUS "Lemon not available, so I/O will be essentially serial.")
  ENDIF(LEMON_FOUND) 
  INCLUDE_DIRECTORIES(BEFORE SYSTEM arch/mpi)
  ADD_SUBDIRECTORY(arch/mpi)
ELSE(USE_MPI) #Scalar build
  INCLUDE_DIRECTORIES(BEFORE SYSTEM arch/scalar)
  ADD_SUBDIRECTORY(arch/scalar)
ENDIF(USE_MPI)

ADD_SUBDIRECTORY(src)
IF(BUILD_TESTS)
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(test)
ENDIF(BUILD_TESTS)
ADD_SUBDIRECTORY(example EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(applications)

ADD_EXECUTABLE(main src/main.cpp)
TARGET_LINK_LIBRARIES(main Debug)

ADD_EXECUTABLE(simon_test src/simon_test.cpp)
TARGET_LINK_LIBRARIES(simon_test SU3 QCD IO Base Core Dirac Weave Contract Input Smear)
