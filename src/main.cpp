#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>

int main(int argc, char **argv)
{
  std::cout << "Dumb test code that should at least compile under all circumstances." << std::endl;
  Core::Field< QCD::Gauge, 8, 8 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.88");
  std::cout << myfield[0][0];
  myfield[0][0].reunitarize();
  std::cout << myfield[0][0];
  return 0;
}

