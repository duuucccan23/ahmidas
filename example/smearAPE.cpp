// $Id$
#include <string>
#include <iomanip>
#include <iostream>

#include <L0/Tool/IO.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Smear/APE.h>

int main(int argc, char **argv)
{
  //This needs some cleaning up, but it is automatically updated by SVN
  std::string Id = "$Id$";
  std::cout << "Executable tag: " << Id.substr(1,Id.length()-2) << std::endl;

  std::cout << "Reading gauge file conf.88 from test directory.\n";
  Core::Field< QCD::Gauge > myfield = Tool::IO::loadILDG<QCD::Gauge >("../../test/conf.88", 8, 8);

  std::cout << "Showing the gauge link in the x-up direction from point (0,0,0,0).\n";
  std::cout << myfield[0][Base::idx_X];

  double strength = 0.5;
  Smear::APE smearer = Smear::APE(strength); //set up a APE smearer with kappa=0.5

  size_t times = 10;
  std::cout << "APE smearing gauge field " << times << " times with strength " << strength << ".\n";

  smearer.smear(myfield, times);

  std::cout << "\nShowing the now smeared gauge link in the x-up direction from point (0,0,0,0).\n";
  std::cout << myfield[0][Base::idx_X];

  return 0;
}
