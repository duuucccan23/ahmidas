SET(Sources_Gamma Gamma/Gamma.cpp)
SET(Sources_Matrix Matrix/Matrix_Matrix.cpp Matrix/Matrix_cout_operator_lshift.cpp
                   Matrix/Matrix_operators.cpp Matrix/Matrix_outer_product.cpp)
ADD_LIBRARY(Dirac ${Sources_Gamma} ${Sources_Matrix})
