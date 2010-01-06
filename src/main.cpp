// $Id$
#include <string>
#include <iomanip>
#include <iostream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Component.h>
#include <L1/Tool.h>

int main(int argc, char **argv)
{
  //This needs some cleaning up, but it is automatically updated by SVN
  std::string Id = "$Id$";
  std::cout << "Executable tag: " << Id.substr(1,Id.length()-2) << std::endl;
  Core::Field< QCD::Spinor, 4, 4 > mysource = Tool::IO::loadScidac<QCD::Spinor, 4, 4 >("../test/source.00.00");
  Core::Field< QCD::Spinor, 4, 4 > myprop = Tool::IO::loadScidac<QCD::Spinor, 4, 4 >("../test/source.00.00");

  std::complex < double> res = Tool::innerProduct(mysource, myprop);
  std::cout << res << std::endl;
  return 0;
}