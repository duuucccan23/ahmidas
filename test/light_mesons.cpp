#include <vector>
#include <complex>
#include <iostream>

#include <L0/Ahmidas.h>
#include <L1/Tool/IO.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Meson.h>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  Dirac::Gamma< 5 > gamma5;

  const size_t L = 4;
  const size_t T = 4;

  std::vector<std::string> propfilesU;
  Core::Propagator *uProp = new Core::Propagator(L, T);

  std::cout << "Testing simplest pion contraction using point sources.\n";
  const std::string filename_base("../../test/source4x4");

  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss <<  ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    propfilesU.push_back(std::string(filename_base).append("_u").append(oss.str()));
    //std::cout << propfilesU[f] << std::endl;
  }

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

  //charged_pion
  Core::Correlator< Dirac::Matrix > C2 = Contract::light_meson_twopoint(uProp, 0, gamma5, gamma5);

  delete uProp;

  std::cout <<  "\nWell tested code gives the following result:\n\n"
            <<  "0  +5.412652273e-01  +0.0000000000e+00\n"
            <<  "1  +1.456410538e-02  +0.0000000000e+00\n"
            <<  "2  +1.637160312e-03  +0.0000000000e+00\n"
            <<  "3  +1.443586407e-02  +0.0000000000e+00\n";

  double tolerance = 1.e-9;

  if  (fabs(C2[0].trace().imag()) < tolerance
    && fabs(C2[1].trace().imag()) < tolerance
    && fabs(C2[2].trace().imag()) < tolerance
    && fabs(C2[3].trace().imag()) < tolerance
    && fabs(C2[0].trace().real()/5.412652273e-01 - 1) < tolerance
    && fabs(C2[1].trace().real()/1.456410538e-02 - 1) < tolerance
    && fabs(C2[2].trace().real()/1.637160312e-03 - 1) < tolerance
    && fabs(C2[3].trace().real()/1.443586407e-02 - 1) < tolerance)
  {
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
