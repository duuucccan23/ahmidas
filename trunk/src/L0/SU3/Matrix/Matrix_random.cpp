#include "Matrix.ih"

SU3::Matrix SU3::Matrix::random()
{
  Base::Ranlux &mything(Base::Ranlux::instance());
  double data[18];
  for (size_t ctr = 0; ctr < 18; ++ctr)
    data[ctr] = mything();
  return Matrix(data);
}
