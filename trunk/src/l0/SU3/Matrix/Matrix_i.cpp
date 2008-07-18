#include "Matrix.ih"

void SU3::Matrix::givens(std::complex< double > &c, std::complex< double > &s, std::complex< double > const &f, std::complex < double > const &g) const
{
  if (g.real() == 0 && g.imag() == 0)
  {
    c = std::complex< double > (1, 0);
    s = std::complex< double > (0, 0);
    return;
  }
  if (f.real() == 0 && f.imag() == 0)
  {
    c = std::complex< double > (0, 0);
    s = sign(conj(g));
    return;
  }
  double inv_sqrt_nrm = 1 / (sqrt(norm(f) + norm(g)));
  c = abs(f) * inv_sqrt_nrm;
  s = sign(f) * conj(g) * inv_sqrt_nrm;
}
