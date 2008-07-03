#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>
#include <QCD/Spinor/Spinor.h>
#include <QCD/Gauge/Gauge.h>
#include <Core/Core.h>
#include <Core/Grid/Grid.h>
#include <Core/Field/Field.h>
#include <Smearing/APE/APE.h>

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
