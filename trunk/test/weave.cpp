#include <iostream>
#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>

int main(int argc, char **argv)
{
  std::cout << "Weaving\n";
  Base::Weave< 4, 4> myWeave = Base::Weave< 4, 4 >::instance();
  return 0;
}
