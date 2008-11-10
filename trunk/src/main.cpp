#include <iostream>
#include <iomanip>

#include <L0/Base/ScidacChecksum.h>

int main(int argc, char **argv)
{
  Base::ScidacChecksum check;
  for (size_t idx = 0; idx < 256; ++idx)
    std::cout << std::dec << idx << "  " << std::hex << check.d_crcTable[idx] << std::endl;
  return 0;
}
