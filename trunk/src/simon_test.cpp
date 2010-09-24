#include <L0/Print.h>
#include <L0/Base/Random.h>
#include <L0/Base/Z2.h>
#include <L0/Ahmidas.h>


int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  size_t samples = 100;
  Base::Z2::instance(1.0, 123456789);
  std::ostringstream ostr("Generating ", std::ios::ate);
  ostr << samples << " random numbers.";
  Print(ostr.str());
  while (samples--)
  {
    std::cout << Base::Random::Z2() <<"\n";
  }
  Print("Done");
  return 0;
}

