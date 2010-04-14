#include <L0/Dirac/Gamma.h>
#include <iostream>

int main(int argc, char **argv)
{
  Dirac::Gamma< -1 > ident = Dirac::Gamma< -1 >();
  Dirac::Gamma< 1 > g1 = Dirac::Gamma< 1 >();
  Dirac::Gamma< 2 > g2 = Dirac::Gamma< 2 >();
  Dirac::Gamma< 3 > g3 = Dirac::Gamma< 3 >();
  Dirac::Gamma< 4 > g4 = Dirac::Gamma< 4 >();
  Dirac::Gamma< 5 > g5 = Dirac::Gamma< 5 >();
  Dirac::Gamma< 15 > g15 = Dirac::Gamma< 15 >();
  Dirac::Gamma< 25 > g25 = Dirac::Gamma< 25 >();
  Dirac::Gamma< 35 > g35 = Dirac::Gamma< 35 >();
  Dirac::Gamma< 45 > g45 = Dirac::Gamma< 45 >();
  Dirac::Gamma< 51 > g51 = Dirac::Gamma< 51 >();
  Dirac::Gamma< 52 > g52 = Dirac::Gamma< 52 >();
  Dirac::Gamma< 53 > g53 = Dirac::Gamma< 53 >();
  Dirac::Gamma< 54 > g54 = Dirac::Gamma< 54 >();
  std::cout << "Gamma< -1 >: identity\n" << ident;
  std::cout << "Gamma< 1 >\n" << g1;
  std::cout << "Gamma< 2 >\n" << g2;
  std::cout << "Gamma< 3 >\n" << g3;
  std::cout << "Gamma< 4 >\n" << g4;
  std::cout << "Gamma< 5 >\n" << g5;
  std::cout << "Gamma< 15 >\n" << g15;
  std::cout << "Gamma< 25 >\n" << g25;
  std::cout << "Gamma< 35 >\n" << g35;
  std::cout << "Gamma< 45 >\n" << g45;
  std::cout << "Gamma< 51 >\n" << g51;
  std::cout << "Gamma< 52 >\n" << g52;
  std::cout << "Gamma< 53 >\n" << g53;
  std::cout << "Gamma< 54 >\n" << g54;
  return 0;
}
