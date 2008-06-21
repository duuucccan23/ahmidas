SET(Sources_Gauge Gauge/Gauge.h)
SET(Sources_Spinor Spinor/Spinor.h)
SET(Sources_WeylSpinor WeylSpinor/WeylSpinor_a.cpp WeylSpinor/WeylSpinor_b.cpp WeylSpinor/WeylSpinor_c.cpp
                       WeylSpinor/WeylSpinor_d.cpp WeylSpinor/WeylSpinor.h)
ADD_LIBRARY(QCD ${Sources_Gauge} ${Sources_Spinor} ${Sources_WeylSpinor})
