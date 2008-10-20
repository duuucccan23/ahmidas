SET(Sources_Tensor     Tensor/Tensor_leftMultiply_a.cpp Tensor/Tensor_rightMultiply_a.cpp
                       Tensor/Tensor_leftMultiply_b.cpp Tensor/Tensor_rightMultiply_b.cpp Tensor/Tensor_Tensor.cpp)
SET(Sources_WeylSpinor WeylSpinor/WeylSpinor_leftMultiply_a.cpp WeylSpinor/WeylSpinor_leftMultiply_b.cpp
                       WeylSpinor/WeylSpinor_leftMultiply_c.cpp WeylSpinor/WeylSpinor_swap.cpp)
ADD_LIBRARY(QCD ${Sources_Tensor} ${Sources_WeylSpinor})
