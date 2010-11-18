/* Tests some of the most basic MPI functionality, to see if it will initialize and finalize properly */
#include <iostream>
#include <mpi.h>

#include <L0/Ahmidas.h>

int main(int argc, char **argv)
{
  bool c0 = (MPI::Is_initialized() || MPI::Is_finalized());
  std::cout << "Before MPI::Init,     MPI::Is_initialized is false and MPI::Is_finalized is false... " << (c0 ? "fail\n" : "pass\n");
  MPI::Init(argc, argv);
  bool c1 = (!MPI::Is_initialized() || MPI::Is_finalized());
  std::cout << " After MPI::Init,     MPI::Is_initialized is true and MPI::Is_finalized is false... " << (c1 ? "fail\n" : "pass\n");
  Ahmidas::Ahmidas my_ahmidas(&argc,&argv);
  MPI::Finalize();
  bool c2 = (!MPI::Is_initialized() || !MPI::Is_finalized());
  std::cout << " After MPI::Finalize, MPI::Is_initialized is true and MPI::Is_finalized is true... " << (c2 ? "fail\n" : "pass\n");
  return (c0 || c1 || c2);
}
