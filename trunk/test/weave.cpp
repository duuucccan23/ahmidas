#include <iostream>
#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>

int main(int argc, char **argv)
{
  std::cout << "Weaving\n";
  Base::Weave myWeave = Base::Weave(4, 8);
  std::cout << myWeave.localVolume() << std::endl;
  std::cout << myWeave.dim(Base::idx_X) << '\t';
  std::cout << myWeave.dim(Base::idx_Y) << '\t';
  std::cout << myWeave.dim(Base::idx_Z) << '\t';
  std::cout << myWeave.dim(Base::idx_T) << std::endl;
  std::cout << myWeave.localSize(Base::idx_X) << '\t';
  std::cout << myWeave.localSize(Base::idx_Y) << '\t';
  std::cout << myWeave.localSize(Base::idx_Z) << '\t';
  std::cout << myWeave.localSize(Base::idx_T) << std::endl;
  return 0;
}
