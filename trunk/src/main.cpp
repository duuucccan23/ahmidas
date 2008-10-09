#include <iostream>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Path.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.88");
  std::cout << myfield[0][0];
  myfield[0][0].reunitarize();
  std::cout << myfield[0][0];
  Core::Field< SU3::Matrix, 8, 8 > staple(myfield, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP);
  return 0;
}

