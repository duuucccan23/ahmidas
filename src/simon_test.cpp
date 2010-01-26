
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Meson.h>


//#define __PRINT__PROPS__
#define __LOAD__PROPS__

#define __REPOSITORY__PROPS__

int main(int argc, char **argv)
{

  const size_t L = 4;
  const size_t T = 4;

  std::cout << "this is a test ..." << std::endl;

  std::vector<std::string> propfilesU;
  std::vector<std::string> propfilesD;

#ifdef __REPOSITORY__PROPS__
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
  const std::string filename_base("../test/propagators/source");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    std::cout << oss.str() << std::endl;
    //propfilesU.push_back(oss.str());
    propfilesU.push_back(oss.str().append("_u"));
    propfilesD.push_back(oss.str().append("_d"));
  }
#endif


  Core::Propagator *uProp = new Core::Propagator(L, T);
  Core::Propagator *dProp = new Core::Propagator(L, T);

#ifdef __LOAD__PROPS__
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
#else
//   uProp->setToRandom();
//   dProp->setToRandom();
#endif

  std::cout << "average difference between u and d propagator: " << uProp->diff(*dProp) << std::endl;
#ifdef __REPOSITORY__PROPS__
  std::cout << "(should be zero up to precision) " << std::endl;
#endif

  Dirac::Gamma<5> gamma5;

#ifdef __PRINT__PROPS__
  for (size_t t=0; t<T; t++)
  {
    std::cout << "timeslice no. " << t << std::endl;
    std::cout << std::endl;

    Core::Propagator::iterator my_iterator = uProp->begin(t);

    int count(0);
    do
    {
      std::cout << "Element no. " << count++ << std::endl;
      std::cout << std::endl;
      std::cout << *(my_iterator);
    }
    while(++my_iterator != uProp->end(t));
  }

  for (size_t t=0; t<T; t++)
  {
    std::cout << "timeslice no. " << t << std::endl;
    std::cout << std::endl;

    Core::Propagator::iterator my_iterator = dProp->begin(t);

    int count(0);
    do
    {
      std::cout << "Element no. " << count++ << std::endl;
      std::cout << std::endl;
      std::cout << *(my_iterator);
    }
    while(++my_iterator != dProp->end(t));
  }
#endif

  Contract::light_meson_twopoint(*uProp, *dProp, gamma5, gamma5);

#ifdef __REPOSITORY__PROPS__
  std::cout <<  "reliable code gives the following result:" << std::endl;
  std::cout <<  "0  +1.616433453    +0" << std::endl;
  std::cout <<  "1  +0.2010924642   +0" << std::endl;
  std::cout <<  "2  +0.04286512674  +0" << std::endl;
  std::cout <<  "3  +0.2105727402   +0" << std::endl;
#endif

  delete uProp;
  delete dProp;

  std::cout << "programm is going to exit normally now\n" << std::endl;


  return EXIT_SUCCESS;
}
