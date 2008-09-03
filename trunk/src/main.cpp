#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{
  cout << "MPI is " << (MPI::Is_initialized() ? "" : "not ") << "initialized." << endl;
  cout << "MPI is " << (MPI::Is_finalized() ? "" : "not ") << "finalized." << endl;
  MPI::Init(argc, argv);
  cout << "MPI is " << (MPI::Is_initialized() ? "" : "not ") << "initialized." << endl;
  cout << "MPI is " << (MPI::Is_finalized() ? "" : "not ") << "finalized." << endl;
  MPI::Finalize();
  cout << "MPI is " << (MPI::Is_initialized() ? "" : "not ") << "initialized." << endl;
  cout << "MPI is " << (MPI::Is_finalized() ? "" : "not ") << "finalized." << endl;
  return 0;
}

