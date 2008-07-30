#include <iostream>
#include <mpi.h>
#include <L0/IO/Lime/Reader.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  IO::Lime::Reader limeTest(argv[1]);
  // Something more intelligent should be done here, probably.
  // Maybe we should construct a special test file for LIME.
  // This particular file is in fact ILDG and needs more care.
  MPI::Finalize();
  return 0;
}
