#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Core/Component.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Path.h>
#include <L1/Tool.h>
#include <L1/Smear/APE.h>

int main(int argc, char **argv)
{
  Core::Field<QCD::Gauge, 12, 4 > field =
      Tool::IO::loadMILC<QCD::Gauge, 12, 4 >("../test/lat.nf8_m02_12x6_b3600_1");

  Tool::fixCoulombGauge(&field);

  return 0;
}
