SET(Sources_WeylSpinor WeylSpinor/WeylSpinor_leftMultiply_a.cpp
                       WeylSpinor/WeylSpinor_leftMultiply_b.cpp
                       WeylSpinor/WeylSpinor_leftMultiply_c.cpp WeylSpinor/WeylSpinor_swap.cpp
                       WeylSpinor.h)
SET(Sources_Tensor Tensor/Tensor_leftMultiply.cpp Tensor/Tensor_rightMultiply.cpp Tensor.h)
ADD_LIBRARY(QCD ${Sources_WeylSpinor} ${Sources_Tensor})
