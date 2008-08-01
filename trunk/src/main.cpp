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
  
  Core::Grid< 8, 8 > grid;
  Core::Field< QCD::Gauge, 8, 8 > field(grid);
  
  IO::ILDG::Generic io("../test/conf.88");
  field.loadDataFromIO(io);

  MPI::Finalize();
  return 0;
}
