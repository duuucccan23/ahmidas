
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Meson.h>
#include <L1/Tool/IO.h>


//#define __GAMMA_TEST__

//#define __PRINT__PROPS__

// #define __REPOSITORY__PROPS_1__

int main(int argc, char **argv)
{

  Dirac::Gamma5 gamma5;
  Dirac::Gamma0 gamma0;
  Dirac::Unity identity;

#ifdef __GAMMA_TEST__
  std::complex< double > tensor_data [144];
  for (int i=0; i<144; i++)
   tensor_data[i] = i;

  QCD::Tensor t1(tensor_data);
  QCD::Tensor t2(tensor_data);

  t1*=gamma0;

  std::ofstream fout1, fout2, fout3;
  fout1.open("bla1");
  fout2.open("bla2");
  fout3.open("bla3");
  fout1 << t1 << std::endl;
  fout2 << t2*gamma0 << std::endl;
  fout3 << gamma0*t2 << std::endl;
  fout1.close();
  fout2.close();
  fout3.close();
  exit(0);
#endif

  const size_t L = 4;
  const size_t T = 4;



  std::cout << "this is a test ..." << std::endl;

  std::vector<std::string> propfilesU;
  std::vector<std::string> propfilesD;

#ifdef __REPOSITORY__PROPS_1__
  const std::string filename_base("../test/source.9999.01");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
    oss << f%4;
    oss << ".inverted";
    oss.flush();
    std::cout << oss.str() << std::endl;
    propfilesU.push_back(oss.str());
    propfilesD.push_back(oss.str());
  }
#else
  const std::string filename_base("../test/source4x4");
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
    propfilesD.push_back(std::string(filename_base).append("_d").append(oss.str()));
    std::cout << propfilesU[f] << std::endl;
    std::cout << propfilesD[f] << std::endl;
  }
#endif


  Core::Propagator *uProp = new Core::Propagator(L, T);
  Core::Propagator *dProp = new Core::Propagator(L, T);

  if (uProp->load(propfilesU, "Scidac"))
  {
    std::cout << "u quark propagator successfully loaded\n" << std::endl;
  }
  else
  {
    std::cout << "error reading u quark  propagator\n" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (dProp->load(propfilesD, "Scidac"))
  {
    std::cout << "d quark propagator successfully loaded\n" << std::endl;
  }
  else
  {
    std::cout << "error reading d quark  propagator\n" << std::endl;
    exit(EXIT_FAILURE);
  }


  std::cout << "average difference between u and d propagator: " << uProp->diff(*dProp) << std::endl;
#ifdef __REPOSITORY__PROPS_1__
  std::cout << "(should be zero up to precision) " << std::endl;
#endif


#ifdef __PRINT__PROPS__

  Core::Propagator::iterator my_iterator = uProp->begin();

  int count(0);
  do
  {
    std::cout << "Element no. " << count++ << std::endl;
    std::cout << std::endl;
    std::cout << *(my_iterator);
  }
  while(++my_iterator != uProp->end());


  Core::Propagator::iterator my_iterator = dProp->begin();

  int count(0);
  do
  {
    std::cout << "Element no. " << count++ << std::endl;
    std::cout << std::endl;
    std::cout << *(my_iterator);
  }
  while(++my_iterator != dProp->end());

#endif

  std::cout <<  "t:" << std::endl;
  Contract::light_meson_twopoint(uProp, 0, gamma5, gamma5);
//   Contract::light_meson_twopoint(uProp, 0, gamma0, gamma0);
//   Contract::light_meson_twopoint(uProp, 0, identity, identity);

//   Contract::light_meson_twopoint(dProp, 0, gamma5, gamma5);
//   Contract::light_meson_twopoint(dProp, 0, gamma0, gamma0);
//   Contract::light_meson_twopoint(dProp, 0, identity, identity);

#ifdef __REPOSITORY__PROPS_1__
  std::cout <<  "reliable code gives the following result:" << std::endl;
  std::cout <<  "0  +1.616433453    +0" << std::endl;
  std::cout <<  "1  +0.2010924642   +0" << std::endl;
  std::cout <<  "2  +0.04286512674  +0" << std::endl;
  std::cout <<  "3  +0.2105727402   +0" << std::endl;
#else
  std::cout <<  "reliable code gives the following result:" << std::endl;
  std::cout <<  " 0  +0.5412652273    +0" << std::endl;
  std::cout <<  " 1  +0.01456410538   +0" << std::endl;
  std::cout <<  " 2  +0.001637160312  +0" << std::endl;
  std::cout <<  " 3  +0.01443586407   +0" << std::endl;
#endif




  delete uProp;
  delete dProp;

  std::cout << "programm is going to exit normally now\n" << std::endl;


  return EXIT_SUCCESS;
}
