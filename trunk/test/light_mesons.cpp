
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include <L1/Tool/IO.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Meson.h>


int main(int argc, char **argv)
{

  std::cout << "test simplest pion contraction using point sources\n" << std::endl;

  Dirac::Gamma5 gamma5;

  const size_t L = 4;
  const size_t T = 4;

  std::vector<std::string> propfilesU;

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

  Core::Propagator *uProp = new Core::Propagator(L, T);

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

  //charged_pion
  Core::Correlator C2 = Contract::light_meson_twopoint(uProp, 0, gamma5, gamma5);

  delete uProp;

  std::cout <<  "\nreliable code gives the following result:\n" << std::endl;
  std::cout <<  " 0  +0.5412652273    +0" << std::endl;
  std::cout <<  " 1  +0.01456410538   +0" << std::endl;
  std::cout <<  " 2  +0.001637160312  +0" << std::endl;
  std::cout <<  " 3  +0.01443586407   +0\n" << std::endl;

  double tolerance = 1.e-10;

  if  (C2[0].trace().imag() == 0.0
    && C2[1].trace().imag() == 0.0
    && C2[2].trace().imag() == 0.0
    && C2[3].trace().imag() == 0.0
    && C2[0].trace().real() <= (0.5412652273   + tolerance) && C2[0].trace().real() >= (0.5412652273   - tolerance)
    && C2[1].trace().real() <= (0.01456410538  + tolerance) && C2[1].trace().real() >= (0.01456410538  - tolerance)
    && C2[2].trace().real() <= (0.001637160312 + tolerance) && C2[2].trace().real() >= (0.001637160312 - tolerance)
    && C2[3].trace().real() <= (0.01443586407  + tolerance) && C2[3].trace().real() >= (0.01443586407  - tolerance))
  {
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
