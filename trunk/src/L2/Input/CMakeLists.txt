SET(Sources_FileReader FileReader.h FileReader/FileReader.cpp)
ADD_LIBRARY(Input ${Sources_FileReader})
SET_TARGET_PROPERTIES(Input PROPERTIES LINKER_LANGUAGE CXX)
