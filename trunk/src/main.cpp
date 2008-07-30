#include <iostream>
#include <mpi.h>
#include <L0/Core/Core.h>
#include <L0/Core/Grid.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/IO/ILDG/Generic.h>
#include <L1/Smearing/Gauge/APE.h>
#include <L1/Smearing/Gauge/Fuzz.h>
#include <L1/Smearing/Spinor/Jacobi.h>

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
/*  Core::Grid< 8,8 > grid;
  if (!grid.rank())
    cout << "MPI and Grid initialized (on node 0).\n";
  Core::Field< QCD::Gauge, 8, 8 > gaugeField(grid);
  if (!grid.rank())
    cout << "Gauge field initialized (on node 0).\n";
*/
  IO::ILDG::Generic io("../test/conf.88");

  double oneval;
  io.read( &oneval, 1);

  cout << oneval << endl;
  MPI::Finalize();
  return 0;
}
