#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  size_t const bufSize = 1024;
  Lime::Reader limeTest(argv[1]);
  std::cout << "Connected to file " << argv[1] << ", with contents of size " << limeTest.size() << ".\n"
            << "Preparing to read." << std::endl;
  double *buffer = new double[bufSize];
  for (size_t ctr = 0; ctr < (limeTest.size() / bufSize); ++ctr)
  {
    limeTest.read(buffer, bufSize);
    if (limeTest.fail())
    {
      std::cout << "Reader went to failed state." << std::endl;
      return 1;
    }
  }
  limeTest.read(buffer, limeTest.size() % bufSize);
    if (limeTest.fail())
    {
      std::cout << "Reader went to failed state after final read." << std::endl;
      return 1;
    }
  std::cout << "Finished reading file to buffer." << std::endl;
  MPI::Finalize();
  return 0;
}
