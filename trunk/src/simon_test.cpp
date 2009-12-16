
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
// #include <L1/Smear/APE.h>

int main(int argc, char **argv)
{

  std::cout << "this it a test ..." << std::endl;

//   const std::string gf_in = "/afs/ifh.de/group/etmc/scratch/poola/dinter/gauge_fields/16x16x16x32/conf.1500";
//   const std::string gf_out =  "/afs/ifh.de/group/nic/scratch/poolb/dinter/tmp/conf.1500";
//   Core::Field< QCD::Gauge, 16, 32 > myfield = Tool::IO::loadILDG<QCD::Gauge, 16, 32 >(gf_in);
//   Tool::IO::saveILDG(myfield, gf_out);
//   // check!

//   Core::Field< QCD::Gauge, 16, 32 > myfield = Tool::IO::loadILDG<QCD::Gauge, 16, 32 >(gf_in);
//   Tool::IO::saveILDG(myfield, gf_out);
//   // check!

  const size_t T = 8;
  const size_t L = 4;
  const std::string src_ildg = "./ahmidas/test/point_src.48";
  Core::Field< QCD::Spinor, L, T > source_field = Tool::IO::loadILDG<QCD::Spinor, L, T >(src_ildg);

  const std::string gf_ildg = "./ahmidas/test/conf.48";
  Core::Field< QCD::GaugeSpinor, L, T > gauge_field = Tool::IO::loadILDG<QCD::Gauge, L, T >(gf_ildg);

  std::cout << source_field[0];
  std::cout << gauge_field[0];

  return 0;
}