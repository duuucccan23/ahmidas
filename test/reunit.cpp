#include <L0/Ahmidas.h>
#include <L0/Debug.h>
#include <L0/Print.h>
#include <fstream>
#include <L0/SU3/Matrix.h>

int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  std::ofstream *out = 0;
  std::ostringstream ostr("", std::ios::ate);
  if (argc > 1)
    out = new std::ofstream(argv[1]);

  bool problem = false;
  for (size_t ctr = 0; ctr < 100000; ++ctr)
  {
    SU3::Matrix mat = SU3::Matrix::random();
    mat.reunitarize();
    SU3::Matrix copy(mat);

    SU3::hcMatrix hc = mat.dagger();
    mat.rightMultiply(hc);
    double error = 0.0;
    for (size_t idx_x = 0; idx_x < 3; ++idx_x)
      for (size_t idx_y = 0; idx_y < 3; ++idx_y)
        error = std::max((idx_x == idx_y) ? std::abs(1.0 - mat(idx_x, idx_y)) : std::abs(mat(idx_x, idx_y)), error);
    if (out)
      (*out) << error;
    if (std::abs(error) > 3e-15)
    {
      ostr.str("At matrix number ");
      ostr << ctr << ":\nStress test failed on unitarity with total deviation " << error << '!';
      Debug(ostr.str());
      problem = true;
    }

    SU3::Matrix copy2(copy);
    copy.reunitarize();

    error = 0.0;
    for (size_t idx_x = 0; idx_x < 3; ++idx_x)
      for (size_t idx_y = 0; idx_y < 3; ++idx_y)
        error += std::abs(copy2(idx_x, idx_y) - copy(idx_x, idx_y));
    if (out)
      (*out) << "\t" << error << std::endl;
    if (std::abs(error) > 3e-15)
    {
      ostr.str("At matrix number ");
      ostr << ctr << ":\nStress test failed on double reunitarization with total deviation " << error << '!';
      Debug(ostr.str());
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

  Print("All matrices properly reunitarized in stress test.");
  return 0;
}
