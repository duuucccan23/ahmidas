#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>
#include <QCD/Gauge/Gauge.h>
#include <Core/Core.h>
#include <Core/Grid/Grid.h>
#include <Core/Field/Field.h>

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
  field.readFromFile("../test/conf.2448");
  if (!grid.rank())
    cout << "Field read and endianness corrected (on node 0).\n";

/*
  short index[4] = {0,0,0,0}; 
  QCD::Gauge myGauge = field.element(index);
  cout << myGauge[0];

  
  QCD::Gauge myGauge = (*(--field.end()));
  cout << myGauge[0];
*/
  MPI::Finalize();
  return 0;
}