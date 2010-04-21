#include "Vector.ih"

SU3::Vector SU3::Vector::random()
{
  SU3::Vector result;
  result.setToRandom();
  return result;
}
