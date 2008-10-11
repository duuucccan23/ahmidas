#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Path.h>
#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>
#include <L1/Smear/APE.h>
#include <L1/Tool.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.88");

  SU3::Matrix mymat(myfield[0][0].dagger());
  std::cout << mymat;
  Core::Field< SU3::Matrix, 8, 8 > result(SU3::Matrix::identity());
  std::cout << result[0];
  Path::step(result, myfield, Base::idx_X, Base::dir_DOWN);
  std::cout << result[0]; 
  Core::Field< SU3::Matrix, 8, 8 > mytest(Path::square(myfield, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP));
  std::cout << mytest[0];
  return 0;
}

