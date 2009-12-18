
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
// #include <L1/Smear/APE.h>

int main(int argc, char **argv)
{

  const size_t L = 4;
  const size_t T = 8;
  
  std::cout << "this is a test ..." << std::endl;

  std::vector<std::string> propfiles;
  
  const std::string filename_base("source");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss << filename_base;
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    std::cout << oss.str() << std::endl;
    propfiles.push_back(oss.str());
  }
  
  Core::Propagator<L, T> *prop = new Core::Propagator<L, T>();
  prop->loadILDG(propfiles);
  

//   const std::string gf_in = "/afs/ifh.de/group/etmc/scratch/poola/dinter/gauge_fields/16x16x16x32/conf.1500";
//   const std::string gf_out =  "/afs/ifh.de/group/nic/scratch/poolb/dinter/tmp/conf.1500";
//   Core::Field< QCD::Gauge, 16, 32 > myfield = Tool::IO::loadILDG<QCD::Gauge, 16, 32 >(gf_in);
//   Tool::IO::saveILDG(myfield, gf_out);
//   // check!

//   Core::Field< QCD::Gauge, 16, 32 > myfield = Tool::IO::loadILDG<QCD::Gauge, 16, 32 >(gf_in);
//   Tool::IO::saveILDG(myfield, gf_out);
//   // check!

//   const size_t T = 8;
//   const size_t L = 4;
//   const std::string src_ildg = "./ahmidas/test/point_src.48";
//   Core::Field< QCD::Spinor, L, T > source_field = Tool::IO::loadILDG<QCD::Spinor, L, T >(src_ildg);
// 
//   const std::string gf_ildg = "./ahmidas/test/conf.48";
//   Core::Field< QCD::GaugeSpinor, L, T > gauge_field = Tool::IO::loadILDG<QCD::Gauge, L, T >(gf_ildg);
// 
//   std::cout << source_field[0];
//   std::cout << gauge_field[0];

  delete prop;

  return 0;
}
