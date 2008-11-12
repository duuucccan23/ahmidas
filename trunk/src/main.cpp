#include <iomanip>
#include <iostream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Path.h>

int main(int argc, char **argv)
{
  Core::Field<QCD::Gauge, 6, 6 > field = Tool::IO::loadMILC<QCD::Gauge, 6, 6 >(std::string("../test/lat.nf8_m02_12x6_b3600_1"));
  Core::Field<SU3::Matrix, 6, 6 > square = Path::square(field, Base::idx_T, Base::dir_UP, Base::idx_X, Base::dir_UP);
  std::cout << square[4] << std::endl;

  return 0;
}
