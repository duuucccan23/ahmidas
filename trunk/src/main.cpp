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
  Core::Field< QCD::Spinor > src1 = Tool::IO::loadScidac<QCD::Spinor >("../test/source4x4_u.00.inverted", 4, 4);

  std::cout << src1[0] << std::endl;
  return 0;
}
