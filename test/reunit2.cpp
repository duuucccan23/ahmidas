#include <iostream>
#include <L0/SU3/Matrix.h>

int main(int, char **)
{
  SU3::Matrix mat = SU3::Matrix::random();
  mat.reunitarize();
  SU3::Matrix copy(mat);
  mat.reunitarize();

  bool c0 = (std::abs(std::norm(mat(0,0) - copy(0,0))) < 1e-14);
  bool c1 = (std::abs(std::norm(mat(0,1) - copy(0,1))) < 1e-14);
  bool c2 = (std::abs(std::norm(mat(0,2) - copy(0,2))) < 1e-14);
  bool c3 = (std::abs(std::norm(mat(1,0) - copy(1,0))) < 1e-14);
  bool c4 = (std::abs(std::norm(mat(1,1) - copy(1,1))) < 1e-14);
  bool c5 = (std::abs(std::norm(mat(1,2) - copy(1,2))) < 1e-14);
  bool c6 = (std::abs(std::norm(mat(2,0) - copy(2,0))) < 1e-14);
  bool c7 = (std::abs(std::norm(mat(2,1) - copy(2,1))) < 1e-14);
  bool c8 = (std::abs(std::norm(mat(2,2) - copy(2,2))) < 1e-14);

  if (c0 && c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8)
  {
    std::cout << "Reunitarization of a previously reunitarized matrix leaves it invariant.\n";
    return 0;
  }
  std::cout << "Reunitarization of a previously reunitarized matrix does not leave it invariant.\n";
  std::cout << "Offending matrix was:\n" << copy;
  std::cout << "Result of reunitarization was:\n" << mat;
  return 1;
}
