#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Tool/IO.h>
#include <L1/Tool/ScidacChecksum.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge > mygauge(4,8);
  Tool::IO::loadILDG(&mygauge, "../test/conf.48");
  return 0;
}
