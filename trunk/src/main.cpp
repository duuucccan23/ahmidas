#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Smear/APE.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.88");
  //Core::Field< QCD::Gauge, 8, 1 > remifield;
  //Base::IO::loadILDG(&remifield, "../test/conf.88_ape_smear.0");

  //std::cout << "According to Remi:\n" << remifield[0][0];

  Smear::APE ape(0.5);
  ape.smear(myfield, 20);
  //std::cout << "According to us:\n" << myfield[0][0];

  return 0;
}
