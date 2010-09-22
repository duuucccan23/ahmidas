#include <L0/Base/Base.h>
#include <L0/Print.h>
#include <L0/Ahmidas.h>

int main(int argc, char **argv)
{
  Ahmidas start(argc, argv);
  Print(&std::cout, "Checking endianness of system...");
  Print(&std::cout, (Base::bigEndian ? "BIG endian" : "little endian"));
  return 0;
}
