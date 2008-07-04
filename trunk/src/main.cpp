#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>
#include <QCD/Spinor/Spinor.h>
#include <QCD/Gauge/Gauge.h>
#include <Core/Core.h>
#include <Core/Grid/Grid.h>
#include <Core/Field/Field.h>
#include <Smearing/APE/APE.h>
#include <Smearing/Fuzz/Fuzz.h>
#include <Smearing/Jacobi/Jacobi.h>

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid< 24, 48 > grid;
  if (!grid.rank())
    cout << "MPI and Grid initialized (on node 0).\n";
  Field< QCD::Spionr, 24, 48 > spinorField;
  if (!grid.rank())
    cout << "Spinor field initialized (on node 0).\n";
  spinorField.readFromFile< float >("../test/prop.2448", "scidac-binary-data");
  if (!grid.rank())
    cout << "Spinor field read and endianness corrected (on node 0).\n";

  MPI::Finalize();
  return 0;
}
