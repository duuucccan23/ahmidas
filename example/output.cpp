#include <L0/Debug.h>
#include <L0/Print.h>
#include <L0/Ahmidas.h>

int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  Print("Hello, I should be on normal cout.");
  Debug("MESSAGE, should be on std::cerr.");
  return 0;
}
