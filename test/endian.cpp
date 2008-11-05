#include <iostream>
#include <L0/Base/Base.h>

int main(int argc, char **argv)
{
  std::cout << "This system is " << (Base::bigEndian ? "BIG endian\n" : "little endian\n");
  return 0;
}
