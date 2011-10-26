#include <L0/Ahmidas.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Identity.h>
#include <L0/Dirac/Matrix.h>
#include <iostream>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  Dirac::Identity ident = Dirac::Identity();
  Dirac::Gamma< 2 > g2 = Dirac::Gamma< 2 >();
  std::cout << "identity\n" << ident;
  std::cout << "Gamma< 2 >\n" << g2;

  const std::complex <double> Id_tmp[16] = { 1., 0., 0., 0.,   
                 0., 1., 0., 0., 
                 0., 0., 1., 0.,
                 0., 0., 0., 1. };

                 
  Dirac::Matrix diracId(Id_tmp);
  
  std::cout << "diracId" << std::endl << diracId << std::endl;
  std::cout << "diracId*g2" << std::endl << diracId*g2 << std::endl;
  std::cout << "g2*diracId" << std::endl << g2*diracId << std::endl;
  

  return 0;
}
