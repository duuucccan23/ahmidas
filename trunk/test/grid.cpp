#include <iostream>
#include <mpi.h>
#include <L0/Base/Grid.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Base::Grid< 4, 4 > *gridp0;
  {
    Base::Grid< 4, 4 > &grid0 = Base::Grid< 4, 4 >::instance();
    gridp0 = &grid0;
    std::cout << "Grid created, adress is: " << gridp0 << std::endl;
  }
  Base::Grid< 4, 4 > &grid1 = Base::Grid< 4, 4 >::instance();
  std::cout << "Another grid created, adress is: " << &grid1 << std::endl;
  bool c0 = &grid1 != gridp0;

  MPI::Finalize();
  return c0;
}
