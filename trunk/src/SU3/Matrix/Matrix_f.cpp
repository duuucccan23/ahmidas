#include "Matrix.ih"

// Based on a complex QR decomposition, 
// implemented following Bindel et al., ACM Trans. on Math. Soft. 28-206 [2002]
// See documentation on unitarization for further information.

void SU3::Matrix::reunitarize()
{
  // Storage for the angles determining the SU3 matrix
  std::complex< double > c[3];
  std::complex< double > s[3];
  // Some storage for intermediate results
  std::complex< double > tilde[4];

  // Perform a series of three successive Givens' rotations
  givens(c[0], s[0], d_data[0], d_data[3]);
  
  tilde[0] = c[0] * d_data[0] + s[0] * d_data[3];
  tilde[1] = c[0] * d_data[1] + s[0] * d_data[4];
  tilde[2] = -std::conj(s[0]) * d_data[1] + std::conj(c[0]) * d_data[4];

  givens(c[1], s[1], tilde[0], d_data[6]);

  tilde[3] = -std::conj(s[1]) * tilde[1] + std::conj(c[1]) * d_data[7];

  givens(c[2], s[2], tilde[2], tilde[3]);

  // Use the resulting angles to build up the unitary core of the matrix.
  d_data[0] = std::conj(c[0] * c[1]);
  d_data[1] = -std::conj(c[0] * s[2]) * s[1] - s[0] * std::conj(c[2]);
  d_data[2] = -std::conj(c[0]) * s[1] * c[2] + s[0] * s[2];
  d_data[3] = std::conj(s[0] * c[1]);
  d_data[4] = -std::conj(s[0] * s[2]) * s[1] + c[0] * std::conj(c[2]);
  d_data[5] = -std::conj(s[0]) * s[1] * c[2] - c[0] * s[2];
  d_data[6] = std::conj(s[1]);
  d_data[7] = c[1] * std::conj(s[2]);
  d_data[8] = c[1] * c[2];
}
