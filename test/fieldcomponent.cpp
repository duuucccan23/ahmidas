#include <iostream>
#include <mpi.h>
#include <Core/Grid/Grid.h>
#include <Core/Field/Field.h>
#include <Core/Component/Component.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid< 4, 8 > grid;
  Core::Field< double, 4, 8 > field(grid);
  Core::Component< double, 4, 8, double > component = field.component< double >(Core::idx_T);
  MPI::Finalize();
  return 0;
}
