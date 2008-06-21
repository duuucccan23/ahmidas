#include "SU3/Matrix/Matrix.h"
#include "SU3/Vector/Vector.h"
#include <iostream>

int main()
{
  SU3::Matrix test = SU3::Matrix::identity();
  std::cout << test;
  return 0;
}

