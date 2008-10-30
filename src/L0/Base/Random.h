#ifndef GUARD_BASE_RANDOM_H
#define GUARD_BASE_RANDOM_H

namespace Base
{
  namespace Random
  {
    double uniform();
    double symmetric();
    double fastUniform();
    double fastSymmetric();
    void setRange(size_t const min, size_t const max);
    size_t range();
    void setZ2Scale(double const scale);
    double Z2();
  }
}

#endif
