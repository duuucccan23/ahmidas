#include <iostream>
#include <mpi.h>
#include <L0/IO/Lime/Reader.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Core.h>
#include <L0/Core/Grid.h>
#include <L0/Core/Field.h>
#include <L1/Smearing/Gauge/APE.h>

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid<8, 8> grid;
  if (!grid.rank())
    cout << "MPI and Grid initialized (on node 0).\n";
  size_t index = 2048;
  cout << grid.rank(index) << endl;
  MPI::Finalize();
  return 0;
}
