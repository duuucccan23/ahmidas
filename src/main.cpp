#include <L0/Print.h>
#include <L0/Ahmidas.h>

int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  Print("Hello");
  return 0;
}
