
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

int main(int argc, char **argv)
{

  const size_t L = 4;
  const size_t T = 4;

  std::cout << "this is a test ..." << std::endl;

  std::vector<std::string> propfiles;

  const std::string filename_base("../test/my_source");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    std::cout << oss.str() << std::endl;
    propfiles.push_back(oss.str());
  }

  Core::Field< QCD::Spinor> *my_field = new Core::Field< QCD::Spinor>(L, T);

  for (int f=0; f<12; f++)
  {
    size_t lattice_site = 11+L*L*L;
    size_t index = f;
    ((*my_field)[lattice_site])(Base::DiracIndex(index/3), Base::ColourIndex(index%3)) = 1.0;
    Tool::IO::saveScidac(*my_field, propfiles[f]);
    ((*my_field)[lattice_site])(Base::DiracIndex(index/3), Base::ColourIndex(index%3)) = 0.0;
  }

  delete my_field;

  Core::Propagator *prop = new Core::Propagator(L, T);
  if (prop->load(propfiles, "Scidac"))
  {
    std::cout << "Propagator structure successfully loaded\n" << std::endl;
  }
  else
  {
    std::cout << "error reading Propagator structure\n" << std::endl;
    exit(EXIT_FAILURE);
  }

  Dirac::Gamma<5> gamma5;
  (*prop) *= gamma5;

  for (size_t t=0; t<T; t++)
  {
    std::cout << "timeslice no. " << t << std::endl;
    std::cout << std::endl;

    Core::Propagator::iterator my_iterator = prop->begin(t);

    int count(0);
    do
    {
      std::cout << "Element no. " << count++ << std::endl;
      std::cout << std::endl;
      std::cout << *(my_iterator);
    }
    while(++my_iterator != prop->end(t));
  }


  delete prop;
  std::cout << "Propagator structure deleted\n" << std::endl;


  return 0;
}
