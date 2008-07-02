#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>
#include <QCD/Gauge/Gauge.h>
#include <Core/Core.h>
#include <Core/Grid/Grid.h>
#include <Core/Field/Field.h>
//#include <Smearing/APE/APE.h>

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
  if (!grid.rank())
    cout << "Field read and endianness corrected (on node 0).\n";
//  Smearing::APE ape(0.5);

  
  MPI::Finalize();
  return 0;
}
