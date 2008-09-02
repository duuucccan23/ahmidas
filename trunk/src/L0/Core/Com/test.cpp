#include <iostream>
#include <mpi.h>
#include "../Com.h"

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Com mycom = Core::Com::setup();
  MPI::Finalize();
  return 0;
}
