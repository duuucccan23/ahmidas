#include <L0/Base/Random.h>

int main(int argc, char **argv)
{
  double cum;
  for (int ctr = 0; ctr < 1000; ++ctr)
    cum += Base::Random::symmetric();
  return 0;
}

