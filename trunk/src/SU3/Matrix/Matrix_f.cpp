#include "Matrix.ih"

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
  double c =   - H2_diag[0] * H2_diag[1] * H2_diag[2]
               - 2 * std::real(H2_off[0] * std::conj(H2_off[1]) * H2_off[2])
               + H2_diag[0] * H2_off_norm[2] + H2_diag[1] * H2_off_norm[1]
               + H2_diag[2] * H2_off_norm[0];

  double Q  = a_3_2 - fac_1_3 * b;
  double sqrt_Q = std::sqrt(Q);
  double sqrt_min_2Q = -2 * sqrt_Q;

  // theta_3 = 1/3 * acos(sqrt(R)/sqrt(Q^3))
  double theta_3 = fac_1_3 * std::acos(std::sqrt(a_3_2 * a_3 + 0.5 * (detH2 - a_3 * b)) / (Q * sqrt_Q)));

  double lambda[3] = { sqrt_min_2Q * std::cos(theta_3) - a_3,
                       sqrt_min_2Q * std::cos(theta_3 + fac_2pi_3) - a_3,
                       sqrt_min_2Q * std::cos(theta_3 - fac_2pi_3) - a_3 };

  // We now know the eigenvalues, time to determine the associated eigenvectors
  // \bar(b) * c - e * a
  std::complex< double > v2_num_zero = std::conj(H2_off[0]) * H2_off[1] - H2_off[2] * H2_diag[0];
  // a*d - b\bar(b)
  std::complex< double > v2_den_zero = H2_diag[0] * H2_diag[1] - H2_off_norm[0];
  // -(a + d)
  std::complex< double > v2_den_one  = -H2_diag[0] - H2_diag[1];

  Matrix eigenVectors;
  for (size_t ctr = 0; ctr < 3; ++ctr)
  {
    size_t offset = 3 * ctr;
    eigenVectors.d_data[offset + 1] = (v2_num_zero + H2_off[2] * lambda[ctr]) / (v2_den_zero + v2_den_one * lambda[ctr] + lambda_2);
    eigenVectors.d_data[offset] = (H2_off[0] * eigenVectors[offset + 1] + H2_off[1]) / (lambda[ctr] - H2_diag[0]);
    double normal = std::pow(1 + eigenVectors[offset] * std::conj(eigenVectors[offset]) + eigenVectors[offset + 1] * std::conj(eigenVectors[offset + 1]), fac_1_3);
    eigenVectors.d_data[offset]     *= normal;
    eigenVectors.d_data[offset + 1] *= normal;
    eigenVectors.d_data[offset + 2]  = normal;
  }

  for (size_t ctr = 0; ctr < 3; ++ctr)
    lambda[ctr] = 1 / std::sqrt(lambda[ctr]);

  // NOTE I'm not sure if it makes sense to store an array of conjugated eigenvectors...
  // In any case: the following will perform the contraction U'*D-(1/2)*U
  Matrix H;
  H.d_data[0] =   lambda[0] * eigenVectors.d_data[0] * std::conj(eigenVectors.d_data[0])
                + lambda[1] * eigenVectors.d_data[3] * std::conj(eigenVectors.d_data[3])
                + lambda[2] * eigenVectors.d_data[6] * std::conj(eigenVectors.d_data[6]);
  H.d_data[1] =   lambda[0] * eigenVectors.d_data[0] * std::conj(eigenVectors.d_data[1])
                + lambda[1] * eigenVectors.d_data[3] * std::conj(eigenVectors.d_data[4])
                + lambda[2] * eigenVectors.d_data[6] * std::conj(eigenVectors.d_data[7]);
  H.d_data[2] =   lambda[0] * eigenVectors.d_data[0] * std::conj(eigenVectors.d_data[2])
                + lambda[1] * eigenVectors.d_data[3] * std::conj(eigenVectors.d_data[5])
                + lambda[2] * eigenVectors.d_data[6] * std::conj(eigenVectors.d_data[8]);
  H.d_data[3] =   lambda[0] * eigenVectors.d_data[1] * std::conj(eigenVectors.d_data[0])
                + lambda[1] * eigenVectors.d_data[4] * std::conj(eigenVectors.d_data[3])
                + lambda[2] * eigenVectors.d_data[7] * std::conj(eigenVectors.d_data[6]);
  H.d_data[4] =   lambda[0] * eigenVectors.d_data[1] * std::conj(eigenVectors.d_data[1])
                + lambda[1] * eigenVectors.d_data[4] * std::conj(eigenVectors.d_data[4])
                + lambda[2] * eigenVectors.d_data[7] * std::conj(eigenVectors.d_data[7]);
  H.d_data[5] =   lambda[0] * eigenVectors.d_data[1] * std::conj(eigenVectors.d_data[2])
                + lambda[1] * eigenVectors.d_data[4] * std::conj(eigenVectors.d_data[5])
                + lambda[2] * eigenVectors.d_data[7] * std::conj(eigenVectors.d_data[8]);
  H.d_data[6] =   lambda[0] * eigenVectors.d_data[2] * std::conj(eigenVectors.d_data[0])
                + lambda[1] * eigenVectors.d_data[5] * std::conj(eigenVectors.d_data[3])
                + lambda[2] * eigenVectors.d_data[8] * std::conj(eigenVectors.d_data[6]);
  H.d_data[7] =   lambda[0] * eigenVectors.d_data[2] * std::conj(eigenVectors.d_data[1])
                + lambda[1] * eigenVectors.d_data[5] * std::conj(eigenVectors.d_data[4])
                + lambda[2] * eigenVectors.d_data[8] * std::conj(eigenVectors.d_data[7]);
  H.d_data[8] =   lambda[0] * eigenVectors.d_data[2] * std::conj(eigenVectors.d_data[2])
                + lambda[1] * eigenVectors.d_data[5] * std::conj(eigenVectors.d_data[5])
                + lambda[2] * eigenVectors.d_data[8] * std::conj(eigenVectors.d_data[8]);

  // We're almost done! All that remains is finishing off the matrix...
  rightMultiply(H);
  operator/=(std::pow(det(d_data), fac_1_3));
}
