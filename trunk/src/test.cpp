#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>
#include <QCD/Gauge/Gauge.h>
#include <Core/Core.h>
#include <Core/Grid/Grid.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Lime::Reader read("../../test/conf.save");
  std::cout << "I managed to connect!" << std::endl;
  std::cout << "Failure state: " << read.fail() << std::endl;
  std::cout << "Size of file: " << read.size() << std::endl;

  QCD::Gauge data[5];
  read.read(data, 5);

  for (size_t ctr = 0; ctr < 5; ++ctr)
    std::cout << data[ctr][0] << std::endl;

  Core::swapEndian(data, data + 5, sizeof(double));

  std::cout << "****  SWAPPED  ****" << std::endl << std::endl;

  for (size_t ctr = 0; ctr < 5; ++ctr)
    std::cout << data[ctr][0] << std::endl;

  std::cout << "****  REUNITARIZED  ****" << std::endl << std::endl;

  for (size_t ctr = 0; ctr < 5; ++ctr)
    data[ctr].reunitarize();
  for (size_t ctr = 0; ctr < 5; ++ctr)
    std::cout << data[ctr][0] << std::endl;
  MPI::Finalize();

  return 0;
}

