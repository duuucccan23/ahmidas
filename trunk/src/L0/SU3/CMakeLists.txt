SET(Sources_Vector Vector/Vector.cpp Vector/Vector_a.cpp Vector/Vector_b.cpp Vector/Vector.h)
SET(Sources_Matrix Matrix/Matrix.cpp Matrix/Matrix_a.cpp Matrix/Matrix_b.cpp Matrix/Matrix_c.cpp
                   Matrix/Matrix_d.cpp Matrix/Matrix_e.cpp Matrix/Matrix_f.cpp Matrix/Matrix_g.cpp
                   Matrix/Matrix_h.cpp Matrix/Matrix_i.cpp Matrix/Matrix_j.cpp Matrix/Matrix_k.cpp
		   Matrix/Matrix.h)
ADD_LIBRARY(SU3 ${Sources_Vector} ${Sources_Matrix})
