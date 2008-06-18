#include "Matrix.ih"

// Based on a complex QR decomposition, 
// implemented following Bindel et al., ACM Trans. on Math. Soft. 28-206 [2002]
// See documentation on unitarization for further information.

void SU3::Matrix::reunitarize()
{
  /* The following code has been commented out, because this procedure contains
     a serious flaw.
  
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
  */
  
  // We calculate the distinct elements of H2, using hermiticity to avoid calculations left and right.
  
  double H2_diag[3] = 
  {
    std::norm(d_data[0]) + std::norm(d_data[3]) + std::norm(d_data[6]),
    std::norm(d_data[1]) + std::norm(d_data[4]) + std::norm(d_data[7]),
    std::norm(d_data[2]) + std::norm(d_data[5]) + std::norm(d_data[8])
  };
  
  std::complex< double > H2_off[3] = 
  {
    d_data[1] * std::conj(d_data[0]) + d_data[4] * std::conj(d_data[3]) + d_data[7] * std::conj(d_data[6]),
    d_data[2] * std::conj(d_data[0]) + d_data[5] * std::conj(d_data[3]) + d_data[8] * std::conj(d_data[6]),
    d_data[2] * std::conj(d_data[1]) + d_data[5] * std::conj(d_data[4]) + d_data[8] * std::conj(d_data[7]),
  };
  
  double H2_off_norm[3] = {std::norm(H2_off[0]), std::norm(H2_off[1]), std::norm(H2_off[2])};  
  
  double detH2 =   H2_diag[0] * H2_diag[1] * H2_diag[2] 
                 + 2 * std::real(H2_off[0] * std::conj(H2_off[1]) * H2_off[2])
                 - H2_diag[0] * H2_off_norm[2] - H2_diag[1] * H2_off_norm[1] - H2_diag[2] * H2_off_norm[0];
  double trH2  =   H2_diag[0] + H2_diag[1] + H2_diag[2];
  double trH2_2 =  trH2 * trH2;
  double c2H2  =   H2_diag[0] * (H2_diag[1] + H2_diag[2]) + H2_diag[1] * H2_diag[2] 
                 - H2_off_norm[0] - H2_off_norm[1] - H2_off_norm[2];
  
  double pfac  = 2 * std::sqrt(trH2_2 - 3 * c2H2); 
  double tabc  = (2 * trH2 * (trH2_2 - 9 * c2H2) + 27 * detH2) / (pfac * (trH2_2 - 3 * c2H2));
  double fac3  = 1.0 / 3.0;
  
  pfac *= fac3;
  
  double r0 = pfac * std::cos(fac3 * std::acos(tabc)) + fac3 * trH2;
  double r1 = -pfac * std::cos(fac3 * std::acos(-tabc)) + fac3 * trH2;
  double r2 = trH2 - r0 - r1;
  
  // TODO Find eigenvectors, etc...
}
