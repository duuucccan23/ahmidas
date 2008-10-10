SET(Sources_Ranlux Ranlux/Ranlux.cpp Ranlux/Ranlux_operator_functional.cpp 
                   Ranlux/Ranlux_initialize.cpp Ranlux/Ranlux_instance.cpp)
SET(Sources_Random Random/Random_symmetric.cpp Random/Random_uniform.cpp)             
ADD_LIBRARY(Base ${Sources_Ranlux} ${Sources_Random})
