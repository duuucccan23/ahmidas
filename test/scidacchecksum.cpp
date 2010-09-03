#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge > mygauge(4,8);
  Tool::IO::loadILDG(&mygauge, "../../test/conf.48");
  std::cout << "In headers: 0x1c231a18 0x3c553b.\n\n";

  Core::Field< QCD::Spinor > myspin(4,4);
  Tool::IO::loadScidac(&myspin, "../../test/source4x4_d.00.inverted");
  std::cout << "In headers: 0xb3c083ec 0xd4211a35.\n";

  std::cout << "This test currently does not get the numbers from the loadScidac function.\n";
  std::cout << "The real version should perform an actual comparison, instead of relying on a user to check.\n";
  return 0;
}
