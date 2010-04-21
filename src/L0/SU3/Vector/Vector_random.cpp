#include "Vector.ih"

SU3::Vector SU3::Vector::random()
{
  GeneralVector< 1 > result;
  result.setToRandom();
  return result;
}
