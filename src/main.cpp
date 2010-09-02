#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  if (argc != 4)
  {
    std::cerr << "This program expects the following input:\n";
    std::cerr << argv[0] << " [L/a] [T/a] [filename]\n";
    std::cerr << "And will then load the file [filename] and safe it with etmc propagator headers.\n";
    exit(EXIT_FAILURE);
  }
  size_t L = atoi(argv[1]);
  size_t T = atoi(argv[2]);
  Core::Field< QCD::Spinor > myprop(L,T);
  Tool::IO::loadScidacUnsafe(&myprop, argv[3], 32);
  Tool::IO::saveScidac(myprop, argv[3]);
  return 0;
}
