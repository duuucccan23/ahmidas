#include "Vector.ih"

  inline SU3::Vector SU3::Vector::random()
{
  Vector result;
  result.setToRandom();
  return result;
}
