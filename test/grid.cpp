#include <iostream>
#include <mpi.h>
#include <l0/IO/Lime/Reader/Reader.h>
#include <l0/QCD/Spinor/Spinor.h>
#include <l0/QCD/Gauge/Gauge.h>
#include <l0/Core/Core.h>
#include <l0/Core/Grid/Grid.h>
#include <l0/Core/Field/Field.h>
#include <l1/Smearing/Gauge/APE/APE.h>

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
