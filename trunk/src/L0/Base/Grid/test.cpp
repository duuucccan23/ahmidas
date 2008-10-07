#include <iostream>
#include <mpi.h>
#include "../Grid.h"

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Base::Grid<4,4> mygrid = Base::Grid<4,4>::instance();
  MPI::Finalize();
  return 0;
}
