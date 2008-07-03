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
  Core::Field< QCD::Gauge, 24, 48 > gaugeField(grid);
  if (!grid.rank())
    cout << "Field initialized (on node 0).\n";
  gaugeField.readFromFile< double >("../test/conf.2448", "ildg-binary-data");
  if (!grid.rank())
    cout << "Field read and endianness corrected (on node 0).\n";
  Core::Field< QCD::Spinor, 24, 48 > spinorField(grid);
  if (!grid.rank())
    cout << "Spinor field initialized (on node 0).\n";
  gaugeField.readFromFile< float >("../test/prop.2448", "scidac-binary-data");
  if (!grid.rank())
    cout << "Spinor field read and endianness corrected (on node 0).\n";
  
  Smearing::APE ape(0.2);
  Smearing::Fuzz fuzz(3);
  Smearing::Jacobi jacobi(0.2);
  
  ape.smear(gaugeField);
  if (!grid.rank())
    cout << "Performed APE smearing.\n";
  fuzz.smear(gaugeField);
  if (!grid.rank())
    cout << "Performed fuzzing.\n";
  jacobi.smear(spinorField, gaugeField);
  if (!grid.rank())
    cout << "Performed Jacobi smearing.\n";

  MPI::Finalize();
  return 0;
}
