
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
#include <L0/QCD/Gauge.h>

int main(int argc, char **argv)
{

  const size_t L = 4;
  const size_t T = 4;

  std::cout << "this is a test ..." << std::endl;

  std::vector<std::string> propfiles;

  const std::string filename_base("../test/source");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
//     oss << f;
    oss << 0 << ".";
    oss.fill('0');
    oss.width(2);
    oss << 0;
//     oss << ".inverted";
    std::cout << oss.str() << std::endl;
    propfiles.push_back(oss.str());
  }

  Core::Propagator *prop = new Core::Propagator(L, T);
  //prop->load(propfiles, "ILDG");
  if (prop->load(propfiles, "Scidac"))
    std::cout << "Propagator structure successfully loaded\n" << std::endl;

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

 // QCD::Tensor::iterator tmp_iterator = my_iterator->begin(Base::col_GREEN,  QCD::ColourStrideSource);

//   std::cout << *(tmp_iterator) << std::endl;
//   while(++tmp_iterator != my_iterator->end(Base::col_GREEN,  QCD::ColourStrideSource))
//     std::cout << *(tmp_iterator) << std::endl;


  delete prop;
  std::cout << "Propagator structure deleted\n" << std::endl;


  return 0;
}
