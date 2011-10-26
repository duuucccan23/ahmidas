#include <L0/Ahmidas.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Identity.h>
#include <L0/Dirac/Matrix.h>
#include <iostream>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  Dirac::Gamma< 25 > g25 = Dirac::Gamma< 25 >();
  std::cout << "Gamma< 25 >\n" << g25;

  const std::complex <double> mdb[16] = { 1., 2., 3., 4.,   
                 5., 6., 7., 8., 
                 9., 10., 11., 12.,
                 13., 14., 15., 16. };

                 
  Dirac::Matrix M(mdb);
  
  std::cout << "M" << std::endl << M << std::endl;
  std::cout << "M*g25" << std::endl << M*g25 << std::endl;
  std::cout << "g25*M" << std::endl << g25*M << std::endl;
  M*=g25;
  std::cout << "M*=g25" << std::endl << M << std::endl;
  Dirac::Matrix M2(mdb);
  M2.left_multiply(g25);
  std::cout << "M.left_multiply(g25)" << std::endl << M2 << std::endl;
  Dirac::Matrix M3(mdb);
  M3.right_multiply(g25);
  std::cout << "M.right_multiply(g25)" << std::endl << M3 << std::endl;
  return 0;
}
