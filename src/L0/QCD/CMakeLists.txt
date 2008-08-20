SET(Sources_Gauge Gauge.h)
SET(Sources_Spinor Spinor.h)
SET(Sources_WeylSpinor WeylSpinor/WeylSpinor_leftMultiply_a.cpp
                       WeylSpinor/WeylSpinor_leftMultiply_b.cpp
                       WeylSpinor/WeylSpinor_leftMultiply_c.cpp WeylSpinor/WeylSpinor_swap.cpp
                       WeylSpinor.h)
ADD_LIBRARY(QCD ${Sources_Gauge} ${Sources_Spinor} ${Sources_WeylSpinor})
