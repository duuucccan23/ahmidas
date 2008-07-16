#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  size_t const bufSize = 1024;
  Lime::Reader limeTest(argv[1], argv[2]);
  size_t doubSize = limeTest.size() / sizeof(double);
  std::cout << "Connected to file " << argv[1] << ", with " << doubSize << " data values of type double.\n"
            << "Preparing to read." << std::endl;
  double *buffer = new double[bufSize];
  for (size_t ctr = 0; ctr < (doubSize / bufSize); ++ctr)
  {
    limeTest.read(buffer, bufSize);
    if (limeTest.fail())
    {
      std::cout << "Reader went to failed state." << std::endl;
      return 1;
    }
  }
  limeTest.read(buffer, doubSize % bufSize);
  std::cout << "Remaining read: " << (doubSize % bufSize) << std::endl;
    if (limeTest.fail())
    {
      std::cout << "Reader went to failed state after final read." << std::endl;
      return 1;
    }
  std::cout << "Finished reading file to buffer." << std::endl;
  MPI::Finalize();
  return 0;
}
