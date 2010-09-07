#include "Matrix.ih"

SU3::Matrix SU3::Matrix::random()
{
  Matrix result;
  result.setToRandom();
  return result;
}
