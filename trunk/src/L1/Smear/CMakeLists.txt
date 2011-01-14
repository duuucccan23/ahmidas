SET(Sources_Smear Fuzz/Fuzz.cpp
                  APE/APE_smear_a.cpp
                  APE/APE_smear_b.cpp
                  APE/APE_smear_c.cpp
                  HYP/HYP_smear_a.cpp
                  HYP/HYP_smear_b.cpp
                  Jacobi/Jacobi_smear_a.cpp
                  Jacobi/Jacobi_smear_b.cpp)
ADD_LIBRARY(Smear ${Sources_Smear})
TARGET_LINK_LIBRARIES(Smear Path Tool)
SET_TARGET_PROPERTIES(Smear PROPERTIES LINKER_LANGUAGE CXX)
