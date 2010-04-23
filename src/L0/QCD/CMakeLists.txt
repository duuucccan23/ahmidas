SET(Sources_Tensor     Tensor/Tensor_leftMultiply_a.cpp Tensor/Tensor_leftMultiply_b.cpp
                       Tensor/Tensor_leftMultiply_c.cpp Tensor/Tensor_rightMultiply_a.cpp
                       Tensor/Tensor_rightMultiply_b.cpp Tensor/Tensor_rightMultiply_c.cpp
                       Tensor/Tensor_operator_assign_a.cpp Tensor/Tensor_operator_assign_b.cpp
                       Tensor/Tensor_getDiracMatrix_a.cpp Tensor/Tensor_getDiracMatrix_b.cpp
                       Tensor/Tensor_getDiracMatrix_c.cpp Tensor/Tensor_getDiracMatrix_d.cpp
                       Tensor/Tensor_getDiracMatrix_alternative.cpp
                       Tensor/Tensor_make_sequential.cpp
                       Tensor/Tensor_make_sequential_d_a.cpp
                       Tensor/Tensor_make_sequential_u_a.cpp
                       Tensor/Tensor_multiplyOuterDiracIndices.cpp
                       Tensor/Tensor_Tensor.cpp)
SET(Sources_WeylSpinor WeylSpinor/WeylSpinor_leftMultiply_a.cpp WeylSpinor/WeylSpinor_leftMultiply_b.cpp
                       WeylSpinor/WeylSpinor_leftMultiply_c.cpp WeylSpinor/WeylSpinor_swap.cpp)
ADD_LIBRARY(QCD ${Sources_Tensor} ${Sources_WeylSpinor})
TARGET_LINK_LIBRARIES(QCD SU3 Dirac)