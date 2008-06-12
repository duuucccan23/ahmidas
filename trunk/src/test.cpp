#include "Core/Buffer/Buffer.h"
#include "Core/Field/Field.h"
#include "Core/Grid/Grid.h"
#include "QCD/Spinor/Spinor.h"

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid< 4, 8 > grid;
  Core::Field< QCD::Spinor, 4, 8 > spinor(grid);
  Core::Buffer< QCD::Spinor > buffer(grid);

  spinor.shift(Core::idx_T, Core::dir_UP);

  MPI::Finalize();
  return 0;
}
