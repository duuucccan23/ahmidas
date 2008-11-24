#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield = Tool::IO::loadILDG<QCD::Gauge, 8, 8 >("../test/conf.88");
  
  std::cout << myfield[0][0];
  return 0;
}
