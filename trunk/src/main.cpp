#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Smear/APE.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield = Tool::IO::loadILDG<QCD::Gauge, 8, 8 >("../test/conf.88");
  
  SU3::Matrix mat1(myfield[0][Base::idx_X]);
  SU3::Matrix mat2(myfield[1][Base::idx_X]);
  SU3::Matrix mat3(myfield[2][Base::idx_X]);
  mat3.leftMultiply(mat2);
  mat3.leftMultiply(mat1);
  std::cout << mat3;
  Core::Field< SU3::Matrix, 8, 8 > mysteps = Path::step(myfield, Base::idx_X, Base::dir_UP, 3);
  mysteps.shift(Base::idx_X, Base::dir_DOWN);
  mysteps.shift(Base::idx_X, Base::dir_DOWN);
  mysteps.shift(Base::idx_X, Base::dir_DOWN);
  std::cout << mysteps[0];
  
  return 0;
}
