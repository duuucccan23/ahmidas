#include "Matrix.ih"

// This reunitarization procedure implements the 'French method' in an explicit way.
// It has been written like this to avoid the single dependency on the LAPACK library
// that was present in previous code. Since this implementation was tailored for the
// specific case of 3x3 hermitian matrices when it comes to the determination of eigenvalues,
// it might also be faster.
// See documentation on unitarization for further information.

void SU3::Matrix::reunitarize()
{
  // We calculate the distinct elements of H2, using hermiticity to avoid calculations left and right.

  static double const fac_1_3 = 1.0 / 3.0;
  static double const fac_2pi_3 = fac_1_3 * 4 * std::acos(0.0);

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

  // a_3 = -(1/3) * tr(H2)
  double a_3 = - fac_1_3 * (H2_diag[0] + H2_diag[1] + H2_diag[2]);
  double a_3_2 = a_3 * a_3;

  // b = c2(H2)
  double b =     H2_diag[0] * (H2_diag[1] + H2_diag[2]) + H2_diag[1] * H2_diag[2]
               - H2_off_norm[0] - H2_off_norm[1] - H2_off_norm[2];
  // c = - det(H2)
  double c = - H2_diag[0] * H2_diag[1] * H2_diag[2]
             - 2 * std::real(H2_off[0] * std::conj(H2_off[1]) * H2_off[2])
             + H2_diag[0] * H2_off_norm[2] + H2_diag[1] * H2_off_norm[1]
             + H2_diag[2] * H2_off_norm[0];
  
  // This algorithm won't work for matrices that aren't full rank!
  if (c == 0)
    return;
  
  
  double Q  = a_3_2 - fac_1_3 * b; 
  double R = a_3_2 * a_3 + 0.5 * (c - a_3 * b);

  // Due to precision issues, unitary matrices will sometimes fail here.
  // In that case, we find that R and Q are very close to 0, 
  // so we can explicitly check here.
  
  
  if (Q <= 0 || R <= 0) // We should already be (almost) unitary
  {
    operator*=(std::pow(det(d_data), -fac_1_3));
    return;
  }
  
  double sqrt_Q = std::sqrt(Q);
  double Q_sqrt_Q = Q * sqrt_Q;
  
  if (R > Q_sqrt_Q) // This implies all sort of nastiness
  {
    operator*=(std::pow(det(d_data), -fac_1_3));
    return;
  }
  
  // theta_3 = 1/3 * acos(arg_theta)
  double theta_3 = fac_1_3 * std::acos(R / (Q_sqrt_Q));

  double sqrt_min_2Q = -2 * sqrt_Q;
  double lambda[3] = { sqrt_min_2Q * std::cos(theta_3) - a_3,
                      sqrt_min_2Q * std::cos(theta_3 + fac_2pi_3) - a_3,
                      sqrt_min_2Q * std::cos(theta_3 - fac_2pi_3) - a_3 };

  // We now know the eigenvalues, time to determine the associated eigenvectors
  // b * e - c * d
  std::complex< double > v2_num_zero = H2_off[0] * H2_off[2] - H2_off[1] * H2_diag[1];
  // a*d - b\bar(b)
  std::complex< double > v2_den_zero = H2_diag[0] * H2_diag[1] - H2_off_norm[0];
  // -(a + d)
  std::complex< double > v2_den_one  = -H2_diag[0] - H2_diag[1];

  Matrix eigenVectors;
  for (size_t ctr = 0; ctr < 3; ++ctr)
  {
    eigenVectors.d_data[ctr] = (v2_num_zero + H2_off[1] * lambda[ctr])
                                  / (v2_den_zero + v2_den_one * lambda[ctr] + lambda[ctr]*lambda[ctr]);
    eigenVectors.d_data[ctr + 3] = (std::conj(H2_off[0]) * eigenVectors.d_data[ctr] + H2_off[2]) /
                                   (lambda[ctr] - H2_diag[1]);
    double normal = std::pow(1 + std::norm(eigenVectors.d_data[ctr]) +
                                 std::norm(eigenVectors.d_data[ctr + 3]), -0.5);
    eigenVectors.d_data[ctr]     *= normal;
    eigenVectors.d_data[ctr + 3] *= normal;
    eigenVectors.d_data[ctr + 6]  = normal;
  }

  for (size_t ctr = 0; ctr < 3; ++ctr)
    lambda[ctr] = 1 / std::sqrt(lambda[ctr]);

  // The following will perform the contraction V*D-(1/2)*V'
  Matrix H;
  H.d_data[0] =   lambda[0] * eigenVectors.d_data[0] * std::conj(eigenVectors.d_data[0])
                + lambda[1] * eigenVectors.d_data[1] * std::conj(eigenVectors.d_data[1])
                + lambda[2] * eigenVectors.d_data[2] * std::conj(eigenVectors.d_data[2]);
  H.d_data[1] =   lambda[0] * eigenVectors.d_data[0] * std::conj(eigenVectors.d_data[3])
                + lambda[1] * eigenVectors.d_data[1] * std::conj(eigenVectors.d_data[4])
                + lambda[2] * eigenVectors.d_data[2] * std::conj(eigenVectors.d_data[5]);
  H.d_data[2] =   lambda[0] * eigenVectors.d_data[0] * std::conj(eigenVectors.d_data[6])
                + lambda[1] * eigenVectors.d_data[1] * std::conj(eigenVectors.d_data[7])
                + lambda[2] * eigenVectors.d_data[2] * std::conj(eigenVectors.d_data[8]);
  H.d_data[3] =   lambda[0] * eigenVectors.d_data[3] * std::conj(eigenVectors.d_data[0])
                + lambda[1] * eigenVectors.d_data[4] * std::conj(eigenVectors.d_data[1])
                + lambda[2] * eigenVectors.d_data[5] * std::conj(eigenVectors.d_data[2]);
  H.d_data[4] =   lambda[0] * eigenVectors.d_data[3] * std::conj(eigenVectors.d_data[3])
                + lambda[1] * eigenVectors.d_data[4] * std::conj(eigenVectors.d_data[4])
                + lambda[2] * eigenVectors.d_data[5] * std::conj(eigenVectors.d_data[5]);
  H.d_data[5] =   lambda[0] * eigenVectors.d_data[3] * std::conj(eigenVectors.d_data[6])
                + lambda[1] * eigenVectors.d_data[4] * std::conj(eigenVectors.d_data[7])
                + lambda[2] * eigenVectors.d_data[5] * std::conj(eigenVectors.d_data[8]);
  H.d_data[6] =   lambda[0] * eigenVectors.d_data[6] * std::conj(eigenVectors.d_data[0])
                + lambda[1] * eigenVectors.d_data[7] * std::conj(eigenVectors.d_data[1])
                + lambda[2] * eigenVectors.d_data[8] * std::conj(eigenVectors.d_data[2]);
  H.d_data[7] =   lambda[0] * eigenVectors.d_data[6] * std::conj(eigenVectors.d_data[3])
                + lambda[1] * eigenVectors.d_data[7] * std::conj(eigenVectors.d_data[4])
                + lambda[2] * eigenVectors.d_data[8] * std::conj(eigenVectors.d_data[5]);
  H.d_data[8] =   lambda[0] * eigenVectors.d_data[6] * std::conj(eigenVectors.d_data[6])
                + lambda[1] * eigenVectors.d_data[7] * std::conj(eigenVectors.d_data[7])
                + lambda[2] * eigenVectors.d_data[8] * std::conj(eigenVectors.d_data[8]);

  // We're almost done! All that remains is finishing off the matrix...
  rightMultiply(H);
  operator*=(std::pow(det(d_data), -fac_1_3));
}
