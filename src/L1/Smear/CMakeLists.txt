SET(Sources_Smear Fuzz.h Fuzz/Fuzz.cpp APE.h APE/APE_smear_a.cpp APE/APE_smear_b.cpp APE/APE_smear_c.cpp
                  Jacobi/Jacobi_smear_a.cpp Jacobi/Jacobi_smear_b.cpp)
ADD_LIBRARY(Smear ${Sources_Smear})
TARGET_LINK_LIBRARIES(Smear Path Tool)
SET_TARGET_PROPERTIES(Smear PROPERTIES LINKER_LANGUAGE CXX)
