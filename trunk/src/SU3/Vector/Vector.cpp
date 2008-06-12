#include "Vector.ih"

namespace
{
  double zero[]    = { 0.0, 0.0,
                       0.0, 0.0,
                       0.0, 0.0 };
  double basis_0[] = { 1.0, 0.0,
                       0.0, 0.0,
                       0.0, 0.0 };
  double basis_1[] = { 0.0, 0.0,
                       1.0, 0.0,
                       0.0, 0.0 };
  double basis_2[] = { 0.0, 0.0,
                       0.0, 0.0,
                       1.0, 0.0 };
}
                           
SU3::Vector const SU3::Vector::s_zero(zero);
SU3::Vector const SU3::Vector::s_basis[3] = { SU3::Vector(basis_0), SU3::Vector(basis_1), SU3::Vector(basis_2) };
