#include "Random.ih"

void Base::Random::setZ2Scale(double const scale)
{
  Base::Z2 &generator(Base::Z2::global());
  generator.setScale(scale);
}
