#include <iostream>
#include <L0/IO/Lime/Reader.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Core.h>
#include <L0/Core/Grid.h>
#include <L0/Core/Field.h>
#include <L1/IO/ILDG/Generic.h>
#include <L1/Smearing/Gauge/APE.h>
#include <L1/Smearing/Gauge/Fuzz.h>
#include <L1/Smearing/Spinor/Jacobi.h>

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid< 24, 48 > grid;
  if (!grid.rank())
    cout << "MPI and Grid initialized (on node 0).\n";
  Core::Field< QCD::Gauge, 24, 48 > gaugeField(grid);
  if (!grid.rank())
    cout << "Gauge field initialized (on node 0).\n";
  IO::ILDG::Generic reader("../test/conf.88");

  MPI::Finalize();
  return 0;
}
