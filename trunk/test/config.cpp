#include <iomanip>
#include <iostream>

#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Tool.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge > myfield(8,8);
  Tool::IO::load(&myfield, "../../test/conf.88", Tool::IO::fileILDG);
  std::cout << "Read conf.88 from test directory.\n";
  double stored = 0.5998194411656625;
  double prec = 1e-14;
  double plaqs = Tool::spatialPlaquette(myfield);
  double plaqds = Tool::spatialDownPlaquette(myfield);
  double plaqt = Tool::temporalPlaquette(myfield);
  double plaqdt = Tool::temporalDownPlaquette(myfield);
  double plaq = 0.5 * (plaqt + plaqs);

  std::cout << std::setprecision(14);
  std::cout << "Spatial plaquette value using UP plaquettes: " << plaqs << std::endl;
  std::cout << "Spatial plaquette value using DOWN plaquettes: " << plaqds << std::endl;
  std::cout << "Temporal plaquette value using UP plaquettes: " << plaqt << std::endl;
  std::cout << "Temporal plaquette value using DOWN plaquettes: " << plaqdt << std::endl;

  std::cout << "Summed plaquette value: "  << plaq << std::endl;
  std::cout << "Stored plaquette value: "  << stored << std::endl;

  bool plaqeq = (fabs(plaq/stored - 1) <= prec);
  std::cout << "This differs " << (plaqeq ? "less" : "more") << " then " << prec << " from the stored value.\n";

  return (!plaqeq);
}
