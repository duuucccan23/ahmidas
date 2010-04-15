SET(Sources_Gamma Gamma/Gamma.cpp)
SET(Sources_Matrix Matrix/Matrix_Matrix.cpp Matrix/Matrix_outer_product.cpp)
ADD_LIBRARY(Dirac ${Sources_Gamma} ${Sources_Matrix})
