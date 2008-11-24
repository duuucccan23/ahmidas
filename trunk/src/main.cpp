#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Spinor, 4, 4 > myfield = Tool::IO::loadScidac<QCD::Spinor, 4, 4 >("../test/AHMSource.0.00");
  
  std::cout << myfield[0][0];
  return 0;
}
