#include <iostream>
#include <SU3/Matrix/Matrix.h>

int main(int argc, char **argv)
{
  SU3::Matrix mat = SU3::Matrix::random();
  mat.reunitarize();
  SU3::hcMatrix hc = mat.dagger();
  mat.rightMultiply(hc);
/*  std::cout << std::norm(mat(0,0)) << std::endl;
  std::cout << std::norm(mat(0,1)) << std::endl;
  std::cout << std::norm(mat(0,2)) << std::endl;
  std::cout << std::norm(mat(1,0)) << std::endl;
  std::cout << std::norm(mat(1,1)) << std::endl;
  std::cout << std::norm(mat(1,2)) << std::endl;
  std::cout << std::norm(mat(2,0)) << std::endl;
  std::cout << std::norm(mat(2,1)) << std::endl;
  std::cout << std::norm(mat(2,2)) << std::endl;*/
  return 0;
}
