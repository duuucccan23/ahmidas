
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

  std::vector<std::string> propfilesU;
  std::vector<std::string> propfilesD;

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


//   Core::Field< QCD::Spinor> *my_field = new Core::Field< QCD::Spinor>(L, T);
// 
//   for (int f=0; f<12; f++)
//   {
//     size_t lattice_site = 11+L*L*L;
//     size_t index = f;
//     ((*my_field)[lattice_site])(Base::DiracIndex(index/3), Base::ColourIndex(index%3)) = 1.0;
//     Tool::IO::saveScidac(*my_field, propfiles[f]);
//     ((*my_field)[lattice_site])(Base::DiracIndex(index/3), Base::ColourIndex(index%3)) = 0.0;
//   }
// 
//   delete my_field;

  Core::Propagator *uProp = new Core::Propagator(L, T);
  if (uProp->load(propfilesU, "Scidac"))
  {
    std::cout << "u quark propagator successfully loaded\n" << std::endl;
  }
  else
  {
    std::cout << "error reading u quark  propagator\n" << std::endl;
    exit(EXIT_FAILURE);
  }

  Core::Propagator *dProp = new Core::Propagator(L, T);
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
  std::cout << "average difference between u and u propagator: " << uProp->diff(*uProp) << std::endl;


  Dirac::Gamma<5> gamma5;
//   (*prop) *= gamma5;
// // 
//   for (size_t t=0; t<T; t++)
//   {
//     std::cout << "timeslice no. " << t << std::endl;
//     std::cout << std::endl;
// 
//     Core::Propagator::iterator my_iterator = uProp->begin(t);
// 
//     int count(0);
//     do
//     {
//       std::cout << "Element no. " << count++ << std::endl;
//       std::cout << std::endl;
//       std::cout << *(my_iterator);
//     }
//     while(++my_iterator != uProp->end(t));
//   }
// 
//   for (size_t t=0; t<T; t++)
//   {
//     std::cout << "timeslice no. " << t << std::endl;
//     std::cout << std::endl;
// 
//     Core::Propagator::iterator my_iterator = dProp->begin(t);
// 
//     int count(0);
//     do
//     {
//       std::cout << "Element no. " << count++ << std::endl;
//       std::cout << std::endl;
//       std::cout << *(my_iterator);
//     }
//     while(++my_iterator != dProp->end(t));
//   }


  Contract::light_meson_twopoint(*uProp, *dProp, gamma5, gamma5);

  delete uProp;
  delete dProp;

  std::cout << "programm is going to exit normally now\n" << std::endl;


  return EXIT_SUCCESS;
}
