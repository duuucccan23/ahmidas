CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Ahmidas CXX)

MESSAGE(STATUS "Host system is " ${CMAKE_HOST_SYSTEM})
MESSAGE(STATUS "Target system name is " ${CMAKE_SYSTEM_NAME})
MESSAGE(STATUS "Target system version is " ${CMAKE_SYSTEM_VERSION})
MESSAGE(STATUS "System processor is " ${CMAKE_SYSTEM_PROCESSOR})

# SET(CMAKE_CXX_COMPILER "icpc")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
SET(CMAKE_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/eigen)

# Option for building the tests
OPTION(BUILD_TESTS "Build all tests packaged with Ahmidas" ON)

INCLUDE_DIRECTORIES(BEFORE SYSTEM src)
SET(CMAKE_CXX_FLAGS "-O3 -Wall -pedantic")

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
#   SET(CMAKE_CXX_COMPILER "mpicxx")
  FIND_PACKAGE(MPI REQUIRED)
  FIND_PACKAGE(LEMON)
  IF(LEMON_FOUND)
    MESSAGE(STATUS "Lemon available, so I/O will be implemented over MPI.")
    INCLUDE_DIRECTORIES(BEFORE SYSTEM ${LEMON_INCLUDE_DIR})
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
MESSAGE(STATUS "Testing turned " ${BUILD_TESTS})

ADD_SUBDIRECTORY(example EXCLUDE_FROM_ALL)
ADD_SUBDIRECTORY(applications)

ADD_EXECUTABLE(main src/main.cpp)
TARGET_LINK_LIBRARIES(main Ahmidas Dirac)


ADD_EXECUTABLE(simon_test src/simon_test.cpp)
IF(USE_MPI)
TARGET_LINK_LIBRARIES(simon_test SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas lemon)
ELSE(USE_MPI)
TARGET_LINK_LIBRARIES(simon_test SU3 IO Base Core Dirac Weave Contract Input Smear Ahmidas)
ENDIF(USE_MPI)

MESSAGE(STATUS "Compiler is " ${CMAKE_CXX_COMPILER})
MESSAGE(STATUS "MPI is " ${USE_MPI})
