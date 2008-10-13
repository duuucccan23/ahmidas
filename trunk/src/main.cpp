#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Smear/APE.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 24, 1 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.2448");
  Smear::APE ape(0.5);
  ape.smear(myfield, 20);

  return 0;
}
