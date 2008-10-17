#include <iostream>
#include <L0/Base/Random.h>

int main(int argc, char **argv)
{
  double t = 0.0;
  for (size_t ctr = 0; ctr < 100000; ++ctr)
    std::cout << Base::Random::fastSymmetric() << std::endl;

  return 0;
}
