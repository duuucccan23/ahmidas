#include <vector>
#include <iostream>

#include <L0/Base/Weave.h>
#include <L0/Core/Propagator.h>
#include <L1/Tool/IO.h>
#include <L2/Contract/Meson.h>

int main(int argc, char **argv)
{
  const size_t L = 4;
  const size_t T = 4;

  Base::Weave weave(L,T);

  std::vector< std::string > propfilesU;

  const std::string filename_base("../../test/source.9999.01");
  for (int f=0; f<4; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    if (weave.isRoot())
      std::cout << oss.str() << std::endl;
    propfilesU.push_back(oss.str());
  }

  Core::StochasticPropagator< 4 > *uProp = new Core::StochasticPropagator< 4 >(L, T);

  if (weave.isRoot())
    std::cout << "\nmemory for Propagator structure allocated\n" << std::endl;

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

  if (weave.isRoot())
    std::cout << "propagator successfully loaded\n" << std::endl;


  std::vector< Core::Correlator > C2 = Contract::light_meson_twopoint_stochastic(*uProp, *uProp);

  delete uProp;

  // Tool::printLightMesonCorrelator(C2, "./lightMesonCorrelator.txt");

  // conventional factors
  (C2[0]) *=  1.0/double(L*L*L);
  (C2[3]) *=  1.0/double(L*L*L);
  (C2[9]) *= -1.0/double(L*L*L);

  if (weave.isRoot())
  {
    std::cout << "Ahmidas results for some gamma-combinations:\n" << std::endl;
    std::cout << "Gamma = gamma5,       Gamma' = gamma5" << std::endl;
    std::cout << C2[0] << std::endl;
    std::cout << "Gamma = gamma0gamma5, Gamma' = gamma0gamma5" << std::endl;
    std::cout << C2[3] << std::endl;
    std::cout << "Gamma = gammaigamma0, Gamma' = gammaigamma0 (summed over i)"<< std::endl;
    std::cout << C2[9] << '\n' << std::endl;
  }

  double tolerance = 1.e-6;

  std::cout << "Reference output from well tested routine:\n\n";
  std::cout << "0  +8.418924e-03    +0.000000e+00\n";
  std::cout << "1  +1.047357e-03    +0.000000e+00\n";
  std::cout << "2  +2.232559e-04    +0.000000e+00\n";
  std::cout << "3  +1.096733e-03    +0.000000e+00\n";

  if  (fabs((C2[0])[0].trace().imag()) < tolerance
    && fabs((C2[0])[1].trace().imag()) < tolerance
    && fabs((C2[0])[2].trace().imag()) < tolerance
    && fabs((C2[0])[3].trace().imag()) < tolerance
    && fabs((C2[0])[0].trace().real()/8.418924e-03 - 1) < tolerance
    && fabs((C2[0])[1].trace().real()/1.047357e-03 - 1) < tolerance
    && fabs((C2[0])[2].trace().real()/2.232559e-04 - 1) < tolerance
    && fabs((C2[0])[3].trace().real()/1.096733e-03 - 1) < tolerance)
  {
    if (weave.isRoot())
      std::cout << "Combination Gamma = gamma5, Gamma' = gamma5 works!\n" << std::endl;
  } else {
    return EXIT_FAILURE;
  }

  std::cout << "0  -1.280935e-04    +0.000000e+00\n";
  std::cout << "1  +7.090511e-06    +0.000000e+00\n";
  std::cout << "2  +2.403232e-05    +0.000000e+00\n";
  std::cout << "3  +4.907771e-05    +0.000000e+00\n";

  if  (fabs((C2[3])[0].trace().imag()) < tolerance
    && fabs((C2[3])[1].trace().imag()) < tolerance
    && fabs((C2[3])[2].trace().imag()) < tolerance
    && fabs((C2[3])[3].trace().imag()) < tolerance
    && fabs((C2[3])[0].trace().real()/-1.280935e-04 - 1) < tolerance
    && fabs((C2[3])[1].trace().real()/7.090511e-06 - 1) < tolerance
    && fabs((C2[3])[2].trace().real()/2.403232e-05 - 1) < tolerance
    && fabs((C2[3])[3].trace().real()/4.907771e-05 - 1) < tolerance)
  {
    if (weave.isRoot())
      std::cout << "Combination Gamma = gamma0gamma5, Gamma' = gamma0gamma5 works!\n" << std::endl;
  } else {
    return EXIT_FAILURE;
  }

  std::cout << "0  +2.811082e-03    +0.000000e+00\n";
  std::cout << "1  +2.908005e-04    +0.000000e+00\n";
  std::cout << "2  +5.055152e-05    +0.000000e+00\n";
  std::cout << "3  +1.540258e-04    +0.000000e+00\n";

  if  (fabs((C2[9])[0].trace().imag()) < tolerance
    && fabs((C2[9])[1].trace().imag()) < tolerance
    && fabs((C2[9])[2].trace().imag()) < tolerance
    && fabs((C2[9])[3].trace().imag()) < tolerance
    && fabs((C2[9])[0].trace().real()/2.811082e-03 - 1) < tolerance
    && fabs((C2[9])[1].trace().real()/2.908005e-04 - 1) < tolerance
    && fabs((C2[9])[2].trace().real()/5.055152e-05 - 1) < tolerance
    && fabs((C2[9])[3].trace().real()/1.540258e-04 - 1) < tolerance)
  {
    if (weave.isRoot())
      std::cout << "Combination Gamma = gammaigamma0, Gamma' = gammaigamma0 (summed over i) works!\n" << std::endl;
  }
  else
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
