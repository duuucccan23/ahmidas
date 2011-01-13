// #include "Matrix.ih"

/***********************************************************
 *   The polar decomposition is implemented
 *   following Fast polar decomposition
 *   N.J. Higham /R.S. Schreiber (1998)
 * *********************************************************/

inline void SU3::Matrix::reunitarize()
{
  double check = 0;
  do
  {
    SU3::Matrix old = *this;
    SU3::Matrix inv = inverse();
    double gamma = sqrt(inv.norm() / norm());
    operator*=(0.5 * gamma);
    operator+=((0.5 / gamma) * inv.dagger());
    old -= *this;
    check = sqrt(old.norm());
  }
  while(check > 1E-15);
  operator/=(std::pow(det(), 1.0/3.0));
}
