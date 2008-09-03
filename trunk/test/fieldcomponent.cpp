#include <iostream>
#include <mpi.h>
#include <L0/Core/Field.h>
#include <L0/Core/Component.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Field< double, 4, 8 > field;
  Core::Component< double, 4, 8, double > component = field.component< double >(Core::idx_T);
  MPI::Finalize();
  return 0;
}
