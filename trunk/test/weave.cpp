#include <iostream>
#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>

int main(int argc, char **argv)
{
  std::cout << "Weaving\n";
  Base::Weave< 4, 4> myWeave = Base::Weave< 4, 4 >::instance();
  std::cout << myWeave.localVolume() << std::endl;
  std::cout << myWeave.dim(Base::idx_X) << std::endl;
  std::cout << myWeave.dim(Base::idx_Y) << std::endl;
  std::cout << myWeave.dim(Base::idx_Z) << std::endl;
  std::cout << myWeave.dim(Base::idx_T) << std::endl;  
  return 0;
}