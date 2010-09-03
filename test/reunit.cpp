#include <iostream>
#include <fstream>
#include <L0/SU3/Matrix.h>

int main(int argc, char **argv)
{
  std::ofstream *out = 0;
  if (argc > 1)
    out = new std::ofstream(argv[1]);

  bool problem = false;
  for (size_t ctr = 0; ctr < 1000; ++ctr)
  {
    SU3::Matrix mat = SU3::Matrix::random();
    mat.reunitarize();
    SU3::hcMatrix hc = mat.dagger();
    mat.rightMultiply(hc);
    double error = 0.0;
    for (size_t idx_x = 0; idx_x < 3; ++idx_x)
      for (size_t idx_y = 0; idx_y < 3; ++idx_y)
	  error = std::max((idx_x == idx_y) ? std::abs(1.0 - mat(idx_x, idx_y)) : std::abs(mat(idx_x, idx_y)), error);
    if (out)
      (*out) << error;
    if (std::abs(error) > 5e-10)
    {
      std::cout << "At matrix number " << ctr << ':' << std::endl;
      std::cout << "Stress test failed on unitarity with total deviation " << error << '!' << std::endl;
      problem = true;
    }
    
    SU3::Matrix copy(mat);
    mat.reunitarize();
    
    error = 0.0;
    for (size_t idx_x = 0; idx_x < 3; ++idx_x)
      for (size_t idx_y = 0; idx_y < 3; ++idx_y)
	  error += std::abs(mat(idx_x, idx_y) - copy(idx_x, idx_y));
    if (out)
      (*out) << "\t" << error << std::endl;
    if (std::abs(error) > 4e-15)
    {
      std::cout << "At matrix number " << ctr << ':' << std::endl;
      std::cout << "Stress test failed on double reunitarization with total deviation " << error << '!' << std::endl;
      problem = true;
    }
  }

  if (out)
  {
    out->close();
    delete out;
  }

  if (problem)
    return 1;
  
  std::cout << "All matrices properly reunitarized in stress test." << std::endl;
  return 0;
}
