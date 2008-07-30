#include <iostream>
#include <mpi.h>
#include <L1/IO/ILDG/Generic.h>

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);

  IO::ILDG::Generic io("../test/conf.2448");

  double oneval;
  io.read( &oneval, 1);

  cout << oneval << endl;
  MPI::Finalize();
  return 0;
}
