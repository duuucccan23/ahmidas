#include <iomanip>
#include <iostream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Base/Random.h>
#include <L0/Core/Field.h>
#include <L0/Core/TMatrix.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Tensor.h>
#include <L1/Tool.h>

int main(int argc, char **argv)
{
  Core::Field<QCD::Gauge, 4, 4 > field = Tool::IO::loadILDG<QCD::Gauge, 4, 4 >(std::string("../test/lat.sample.l4444.ildg"));
  Tool::reunitarize(&field);
  std::cerr << field[0][0] << std::endl;

  Core::Field<QCD::Gauge, 4, 4 > sfield = Tool::IO::loadMILC<QCD::Gauge, 4, 4 >(std::string("../test/lat.sample.l4444"));
  Tool::reunitarize(&sfield);
  std::cerr << sfield[0][0] << std::endl;
  std::cerr << det(sfield[0][0]) << std::endl;

  Core::Field< QCD::Gauge, 8, 8 > qField = Tool::IO::loadILDG< QCD::Gauge, 8, 8 >(std::string("../test/conf.88"));
  Tool::reunitarize(&field);
  std::cerr << qField[0][0] << std::endl;

  return 0;
}
