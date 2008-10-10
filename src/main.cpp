#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Path.h>
#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.88");
  Core::Field< SU3::Matrix, 8, 8 > mystaple = Path::staple(myfield, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP);
  std::cout << myfield[1][Base::idx_Y];
  myfield.shift(Base::idx_X, Base::dir_DOWN);
  std::cout << myfield[0][Base::idx_Y];
  return 0;
}

