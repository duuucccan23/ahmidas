#include <iostream>
#include <L0/Ahmidas.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  if (argc != 5)
  {
    std::cerr << "This program expects the following input:\n";
    std::cerr << argv[0] << " [L/a] [T/a] [filename] [32|64]\n";
    std::cerr << "And will then load the file [filename] and save it with ETMC propagator Lime information.\n";
    exit(EXIT_FAILURE);
  }
  size_t L = atoi(argv[1]);
  size_t T = atoi(argv[2]);
  size_t prec = atoi(argv[4]);
  Core::Field< QCD::Spinor > myprop(L,T);
  Tool::IO::loadScidacUnsafe(&myprop, argv[3], prec);
  Tool::IO::saveScidac(myprop, argv[3]);
  return 0;
}
