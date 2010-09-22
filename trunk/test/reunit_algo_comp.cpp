#include <iostream>
#include <limits>
#include <L0/SU3/Matrix.h>
#include <Eigen/Dense>

namespace Eigen
{
  typedef Matrix< std::complex< double >, 3, 3 > Matrix3c;
}

Eigen::Matrix3c toEigen(SU3::Matrix const &mat)
{
  Eigen::Matrix3c mn;
  mn(0, 0) = mat(0, 0);
  mn(1, 0) = mat(1, 0);
  mn(2, 0) = mat(2, 0);
  mn(0, 1) = mat(0, 1);
  mn(1, 1) = mat(1, 1);
  mn(2, 1) = mat(2, 1);
  mn(0, 2) = mat(0, 2);
  mn(1, 2) = mat(1, 2);
  mn(2, 2) = mat(2, 2);

  return mn;
}

SU3::Matrix toSU3(Eigen::Matrix3c const &mat)
{
  SU3::Matrix mn;
  mn(0, 0) = mat(0, 0);
  mn(1, 0) = mat(1, 0);
  mn(2, 0) = mat(2, 0);
  mn(0, 1) = mat(0, 1);
  mn(1, 1) = mat(1, 1);
  mn(2, 1) = mat(2, 1);
  mn(0, 2) = mat(0, 2);
  mn(1, 2) = mat(1, 2);
  mn(2, 2) = mat(2, 2);

  return mn;
}

void reunitarize_new_method_through_eigen(SU3::Matrix &mat_SU3)
{
  Eigen::Matrix3c mat = toEigen(mat_SU3);
  Eigen::Matrix3c mold;
  do
  {
    mold = mat;
    Eigen::Matrix3c minv = mat.inverse();
    std::complex< double > gamma = sqrt(minv.norm() / mat.norm());
    mat = 0.5 * (gamma * mat + minv.conjugate().transpose() / gamma);
  }
  while ((mat - mold).norm() > 1E-15);
  mat /= std::pow(mat.determinant(), 1.0/3.0);
  mat_SU3 = toSU3(mat);
}

void reunitarize_old_method_through_eigen(SU3::Matrix &mat_SU3)
{
  Eigen::Matrix3c mat = toEigen(mat_SU3);
  Eigen::Matrix3c h2 = mat * mat.transpose().conjugate();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3c> es(h2);
  Eigen::Vector3d eigs = es.eigenvalues();
  eigs.array() = eigs.array().sqrt().inverse();
  mat = es.eigenvectors() * eigs.asDiagonal() * es.eigenvectors().transpose().conjugate() * mat;
  mat /= std::pow(mat.determinant(), 1.0/3.0);
  mat_SU3 = toSU3(mat);
}

int main(int argc, char **argv)
{
  double const prec = 1E-16;

  bool problem = false;
  for (size_t ctr = 0; ctr < 100000; ++ctr)
  {
    SU3::Matrix new_method = SU3::Matrix::random();
    SU3::Matrix old_method(new_method);
//     SU3::Matrix new_method_eigen(new_method);
//     SU3::Matrix old_method_eigen(new_method);

    new_method.reunitarize();
    old_method.DEPRECATED_reunitarize();
//     reunitarize_old_method_through_eigen(old_method_eigen);
//     reunitarize_new_method_through_eigen(new_method_eigen);

//     std::cout << "Old method (deprecated):\n" << old_method;
//     std::cout << "New method:\n" << new_method;
//    std::cout << "Eigen version old method:\n" << old_method_eigen;
//    std::cout << "Eigen version new method:\n" << new_method_eigen;
    new_method -= old_method;
//     std::cout << "Difference:\n" << new_method;
//     std::cout << "Norm:\n" << new_method.norm() << "\n\n" << std::endl;

    problem = problem || (new_method.norm() > prec);

  }

  if (problem)
  {
    std::cout << "Problem found: difference in results between algorithms is LARGER than " << prec << '.' << std::endl;
    return 1;
  }

  std::cout << "Success: Difference in results between algorithms is smaller than " << prec << '.' << std::endl;
  return 0;
}
