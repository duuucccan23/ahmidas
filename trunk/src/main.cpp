// $Id$
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
  Core::Field< QCD::Gauge, 20, 48 > myfield = Tool::IO::loadILDG<QCD::Gauge, 20, 48 >("../test/conf.2448");
  Tool::IO::saveILDG(myfield, "../test/conf.2048");
  
  return 0;
}
