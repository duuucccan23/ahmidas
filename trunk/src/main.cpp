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
  Core::Grid<24, 48> grid;
  if (!grid.rank())
    cout << "MPI and Grid initialized (on node 0).\n";
  Core::Field<QCD::Gauge, 24, 48> field(grid);
  if (!grid.rank())
    cout << "Field initialized (on node 0).\n";
  field.readFromFile< double >("../test/conf.2448", "ildg-binary-data");
  cout << grid.dim(0) << " " << grid.dim(1) << " " << grid.dim(2) << " " << grid.dim(3) << endl;
  cout << grid.dimSize(0) << " " << grid.dimSize(1) << " " << grid.dimSize(2) << " " << grid.dimSize(3) << endl;
  if (!grid.rank())
    cout << "Field read and endianness corrected (on node 0).\n";
  Smearing::APE ape(0.5);
  ape.smear(field, 1);

  MPI::Finalize();
  return 0;
}
