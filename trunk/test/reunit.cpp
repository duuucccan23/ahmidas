#include <iostream>
#include <fstream>
#include <L0/SU3/Matrix.h>

int main(int, char **)
{
  srand(time(0));
  
  for (size_t ctr = 0; ctr < 10000; ++ctr)
  {
    SU3::Matrix mat = SU3::Matrix::random(rand());
    mat.reunitarize();
    SU3::hcMatrix hc = mat.dagger();
    mat.rightMultiply(hc);
    double error = -3.0;
    for (size_t idx_x = 0; idx_x < 3; ++idx_x)
      for (size_t idx_y = 0; idx_y < 3; ++idx_y)
	error += std::abs(mat(idx_x, idx_y));
    if (std::abs(error) > 1e-10)
    {
      std::cout << "Stress test failed on unitarity with total deviation " << error << '!' << std::endl;
      return 1;
    }
    
    SU3::Matrix copy(mat);
    mat.reunitarize();
    
    error = 0.0;
    for (size_t idx_x = 0; idx_x < 3; ++idx_x)
      for (size_t idx_y = 0; idx_y < 3; ++idx_y)
	error += std::abs(mat(idx_x, idx_y) - copy(idx_x, idx_y));
    if (std::abs(error) > 1e-10)
    {
      std::cout << "Stress test failed on double reunitarization with total deviation " << error << '!' << std::endl;
      return 1;
    }
  }

  std::cout << "All matrices properly reunitarized in stress test." << std::endl;
  return 0;
}
