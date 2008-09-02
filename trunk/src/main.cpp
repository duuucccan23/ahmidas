#include <iostream>
#include <mpi.h>
#include <L0/Core/Field.h>
#include <L0/Core/Grid.h>
#include <L0/QCD/Gauge.h>
#include <L1/IO/ILDG/Generic.h>

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);

  IO::Lime::Reader io("../test/conf.88");

  cout << "here!\n";

  MPI::Finalize();
  return 0;
}

