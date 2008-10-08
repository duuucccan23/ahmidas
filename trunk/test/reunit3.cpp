#include <iostream>
#include <L0/SU3/Matrix.h>

int main(int, char **)
{
  double vals[18] = { 6.5488492837371E-02,  1.9773604830435E-01,
				    -3.6449121584861E-01, -4.2627838845700E-01,
					 4.3265249921886E-01,  6.7443043861488E-01,
					-2.0268039838887E-01, -6.2614004491101E-01,
					-4.7460514215932E-01, -4.9563949900024E-01,
					-1.9029516948811E-01, -2.4443505232266E-01,
					-1.6377172629791E-01, -7.0474096974707E-01,
					 2.9872624674871E-01,  3.5395898727658E-01,
					 2.5995936491859E-01,  4.4092604895385E-01};
					
  SU3::Matrix mat(vals);
  SU3::Matrix copy(mat);
  mat.reunitarize();

  bool c0 = (std::abs(std::norm(mat(0,0) - copy(0,0))) < 1e-14);
  bool c1 = (std::abs(std::norm(mat(0,1) - copy(0,1))) < 1e-14);
  bool c2 = (std::abs(std::norm(mat(0,2) - copy(0,2))) < 1e-14);
  bool c3 = (std::abs(std::norm(mat(1,0) - copy(1,0))) < 1e-14);
  bool c4 = (std::abs(std::norm(mat(1,1) - copy(1,1))) < 1e-14);
  bool c5 = (std::abs(std::norm(mat(1,2) - copy(1,2))) < 1e-14);
  bool c6 = (std::abs(std::norm(mat(2,0) - copy(2,0))) < 1e-14);
  bool c7 = (std::abs(std::norm(mat(2,1) - copy(2,1))) < 1e-14);
  bool c8 = (std::abs(std::norm(mat(2,2) - copy(2,2))) < 1e-14);

  if (c0 && c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8)
  {
    std::cout << "Reunitarization of known problematic unitary matrix leaves it invariant.\n";
    return 0;
  }
  std::cout << "Reunitarization of a known problematic matrix does not leave it invariant.\n";
  std::cout << "Offending matrix was:\n" << copy;
  std::cout << "Result of reunitarization was:\n" << mat;
  return 1;
}
