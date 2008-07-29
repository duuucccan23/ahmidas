#include <iostream>
#include <L0/IO/Lime/Reader.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Core.h>
#include <L0/Core/Grid.h>
#include <L0/Core/Field.h>
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
  Field< QCD::Spinor, 24, 48 > spinorField;
  if (!grid.rank())
    cout << "Spinor field initialized (on node 0).\n";
  spinorField.readFromFile< float >("../test/prop.2448", "scidac-binary-data");
  if (!grid.rank())
    cout << "Spinor field read and endianness corrected (on node 0).\n";

  MPI::Finalize();
  return 0;
}
