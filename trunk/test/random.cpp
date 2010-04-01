#include <iostream>
#include <L0/Base/Random.h>

int main(int argc, char **argv)
{
  double cum;
  size_t samples = 100000;
  std::cout << "Generating " << samples << " random numbers.\n";
  while (samples--)
    cum += Base::Random::symmetric();
  return 0;
}

