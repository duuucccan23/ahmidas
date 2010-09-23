#include <L0/Print.h>
#include <L0/Base/Random.h>
#include <L0/Ahmidas.h>

int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  double cum;
  size_t samples = 100000;
  std::ostringstream ostr("Generating ", std::ios::ate);
  ostr << samples << " random numbers.";
  Print(ostr.str());
  while (samples--)
    cum = Base::Random::symmetric();
  Print("Done");
  return 0;
}

