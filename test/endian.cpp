#include <iostream>
#include <l0/Core/Core.h>

int main(int argc, char **argv)
{
  std::cout << "This system is " << (Core::big_endian() ? "BIG endian\n" : "little endian\n");
  return 0;
}
