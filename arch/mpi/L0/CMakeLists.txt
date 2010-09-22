ADD_SUBDIRECTORY(Base)

SET(Sources_Print Print/Print_Print.cpp)
ADD_LIBRARY(Print ${Sources_Print})

SET(Sources_Debug Debug/Debug_Debug.cpp)
ADD_LIBRARY(Debug ${Sources_Debug})

SET(Sources_Ahmidas Ahmidas/Ahmidas_Ahmidas.cpp Ahmidas/Ahmidas_destroy.cpp)
ADD_Library(Ahmidas ${Sources_Ahmidas})
