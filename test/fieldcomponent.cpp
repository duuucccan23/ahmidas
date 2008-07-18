#include <iostream>
#include <mpi.h>
#include <l0/Core/Grid/Grid.h>
#include <l0/Core/Field/Field.h>
#include <l0/Core/Component/Component.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid< 4, 8 > grid;
  Core::Field< double, 4, 8 > field(grid);
  Core::Component< double, 4, 8, double > component = field.component< double >(Core::idx_T);
  MPI::Finalize();
  return 0;
}
