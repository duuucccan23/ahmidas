# - Try to find Lemon
# Once done, this will define
#
#  LEMON_FOUND - system has Lemon
#  LEMON_INCLUDE_PATH - the Lemon include directories
#  LEMON_LIBRARY - link these to use Lemon

# Find the include directory
FIND_PATH(LEMON_INCLUDE_PATH NAMES lemon.h PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include $ENV{LEMONDIR}/include PATH_SUFFIXES lemon)

# Find the library itself
FIND_LIBRARY(LEMON_LIBRARY NAMES liblemon.a PATHS /usr/lib /usr/local/lib ${CMAKE_LIBRARY_PATH} ${CMAKE_PREFIX_PATH}/lib $ENV{LEMONDIR}/lib $ENV{LEMONDIR}/src PATH_SUFFIXES lemon)

# Set the include dir variables and the libraries and let libfind_process do the rest.
IF (LEMON_INCLUDE_PATH AND LEMON_LIBRARY)
   SET(LEMON_FOUND TRUE)
ENDIF (LEMON_INCLUDE_PATH AND LEMON_LIBRARY)

IF (LEMON_FOUND)
   IF (NOT LEMON_FIND_QUIETLY)
      MESSAGE(STATUS "Found Lemon: ${LEMON_LIBRARY}")
   ENDIF (NOT LEMON_FIND_QUIETLY)
ELSE (LEMON_FOUND)
   IF (LEMON_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find required library Lemon")
   ENDIF (LEMON_FIND_REQUIRED)
ENDIF (LEMON_FOUND)
