SET(Sources_Ranlux Ranlux/Ranlux.cpp Ranlux/Ranlux_operator_functional.cpp 
                   Ranlux/Ranlux_initialize.cpp Ranlux/Ranlux_instance.cpp)
SET(Sources_Knuth  Knuth/Knuth.cpp Knuth/Knuth_operator_functional.cpp
                   Knuth/Knuth_initialize.cpp Knuth/Knuth_instance.cpp)
SET(Sources_Random Random/Random_symmetric.cpp Random/Random_uniform.cpp
                   Random/Random_fastUniform.cpp Random/Random_fastSymmetric.cpp
                   Random/Random_Z2.cpp Random/Random_setZ2Scale.cpp)
SET(Sources_Z2     Z2/Z2.cpp Z2/Z2_Z2.cpp Z2/Z2_instance.cpp)
ADD_LIBRARY(Base ${Sources_Ranlux} ${Sources_Knuth} ${Sources_Random} ${Sources_Z2})
