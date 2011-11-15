#include <iostream>
#include <cstdlib>

#include <L0/Ahmidas.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Tool.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);
  if (argc != 5) {
    std::cout << "Usage: " << argv[0] << " [L] [T] [infile] [outfile]\n\n";
    return 1;
  }
  int L = atoi(argv[1]);
  int T = atoi(argv[2]);
  if ((L <= 0) | (T <= 0)) {
    std::cout << "Non-positive or non-numeric argument given for L and/or T.\n";
    std::cout << "Unable to proceed.\n";
    std::cout << "Usage: " << argv[0] << " [L] [T] [infile] [outfile]\n\n";
    return 1;
  }    
  Core::Field< QCD::Gauge > myfield(L,T);
  Tool::IO::load(&myfield, argv[3], Tool::IO::fileMILC);
  Tool::IO::saveILDG(myfield, argv[4]);
  return 0;
}
