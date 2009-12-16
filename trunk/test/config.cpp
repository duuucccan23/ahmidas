#include <string>
#include <iomanip>
#include <iostream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Tool.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield = Tool::IO::loadILDG<QCD::Gauge, 8, 8 >("../../test/conf.88");
  std::cout << "Read conf.88 from test directory.\n";
  double stored = 0.5998194411656625;
  double prec = 1e-14;
  double plaqs = Tool::spatialPlaquette(myfield);
  double plaqt = Tool::temporalPlaquette(myfield);
  double plaq = 0.5 * (plaqt + plaqs);

  std::cout << "Summmed plaquette value: " << std::setprecision(14) << plaq << std::endl;
  bool plaqeq = (abs(plaq/stored - 1) <= prec);
  std::cout << "This differs " << (plaqeq ? "less" : "more") << " then " << prec << " from the stored value.\n";

  return (!plaqeq);
}
