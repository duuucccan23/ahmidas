#include <L0/Debug.h>
#include <L0/Print.h>
#include <L0/Ahmidas.h>

int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  Print("Hello, I should be on normal cout.");
  Debug("THIS IS A NORMAL DEBUG TEST MESSAGE");
  return 0;
}
