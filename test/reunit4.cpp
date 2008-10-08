#include <iostream>
#include <L0/SU3/Matrix.h>

int main(int, char **)
{
  srand(time(0));
      
  for (size_t ctr = 0; ctr < 1000; ++ctr)
  {
    SU3::Matrix mat = SU3::Matrix::random(rand());
    SU3::Matrix copy(mat);
    mat.reunitarize();
    SU3::hcMatrix hc = mat.dagger();
    mat.rightMultiply(hc);

    bool c0 = std::abs(std::abs(mat(0,0)) - 1) >= 1e-11;
    bool c1 = std::abs(std::abs(mat(0,1)))     >= 1e-11;
    bool c2 = std::abs(std::abs(mat(0,2)))     >= 1e-11;
    bool c3 = std::abs(std::abs(mat(1,0)))     >= 1e-11;
    bool c4 = std::abs(std::abs(mat(1,1)) - 1) >= 1e-11;
    bool c5 = std::abs(std::abs(mat(1,2)))     >= 1e-11;
    bool c6 = std::abs(std::abs(mat(2,0)))     >= 1e-11;
    bool c7 = std::abs(std::abs(mat(2,1)))     >= 1e-11;
    bool c8 = std::abs(std::abs(mat(2,2)) - 1) >= 1e-11;
    
    if (c0 || c1 || c2 || c3 || c4 || c5 || c6 || c7 || c8)
    {
      std::cout << "Non-SU3 matrix produced during stress test of reunitarization.\n";
      std::cout << "Offending matrix was:\n" << copy;
      std::cout << "Result of reunitarization was:\n" << mat;
      return 1;
    }
  }
  std::cout << "Reunitarization of a random complex 3x3 matrix produces an SU3 matrix.\n";
  return 0;
}
