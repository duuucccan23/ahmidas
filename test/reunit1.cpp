#include <iostream>
#include <l0/SU3/Matrix/Matrix.h>

int main(int, char **)
{
  SU3::Matrix mat = SU3::Matrix::random();
  mat.reunitarize();
  SU3::hcMatrix hc = mat.dagger();
  mat.rightMultiply(hc);

  bool c0 = (std::abs(std::norm(mat(0,0)) - 1) < 1e-12); ///NOTE: maybe an input tolerance parameter can be used here
  bool c1 = (std::abs(std::norm(mat(0,1)) - 0) < 1e-24);
  bool c2 = (std::abs(std::norm(mat(0,2)) - 0) < 1e-24);
  bool c3 = (std::abs(std::norm(mat(1,0)) - 0) < 1e-24);
  bool c4 = (std::abs(std::norm(mat(1,1)) - 1) < 1e-12);
  bool c5 = (std::abs(std::norm(mat(1,2)) - 0) < 1e-24);
  bool c6 = (std::abs(std::norm(mat(2,0)) - 0) < 1e-24);
  bool c7 = (std::abs(std::norm(mat(2,1)) - 0) < 1e-24);
  bool c8 = (std::abs(std::norm(mat(2,2)) - 1) < 1e-12);

  if (c0 && c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8)
  {
    std::cout << "Reunitarization of a random complex 3x3 matrix produces an SU3 matrix.\n";
    return 0;
  }
  std::cout << "Non-SU3 matrix produced when unitarizing a random complex 3x3 matrix.\n";
  return 1;
}
