#include <iostream>
#include <mpi.h>
#include <L0/Core/Grid.h>
#include <L1/QCD/Tensor.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid< 8, 8 > grid;
  QCD::Tensor< 8, 8 > tensor(grid);
  MPI::Finalize();
  return 0;
}
