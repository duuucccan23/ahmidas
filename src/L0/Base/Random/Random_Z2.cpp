#include "Random.ih"

double Base::Random::Z2()
{
  static Base::Z2 &generator(Base::Z2::global());
  return generator();
}
