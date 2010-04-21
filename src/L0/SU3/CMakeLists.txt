SET(Sources_Matrix Matrix/Matrix_givens.cpp Matrix/Matrix_cout_operator_lshift.cpp
                   Matrix/Matrix_leftMultiply_a.cpp Matrix/Matrix_leftMultiply_b.cpp
                   Matrix/Matrix_Matrix.cpp Matrix/Matrix_random.cpp
                   Matrix/Matrix_reunitarize.cpp Matrix/Matrix_rightMultiply_a.cpp
                   Matrix/Matrix_rightMultiply_b.cpp Matrix/Matrix.cpp
                   Matrix/Matrix_operator_addeq.cpp Matrix/Matrix_operator_subeq.cpp
                   Matrix/Matrix_operator_assign.cpp)
SET(Sources_Vector Vector/Vector_leftMultiply_a.cpp Vector/Vector_leftMultiply_b.cpp
                   Vector/Vector.cpp Vector/Vector_cout_operator_lshift.cpp Vector/Vector_random.cpp)
ADD_LIBRARY(SU3 ${Sources_Vector} ${Sources_Matrix})
TARGET_LINK_LIBRARIES(SU3 Base)
