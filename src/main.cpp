// $Id$
#include <string>
#include <iomanip>
#include <iostream>

#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Component.h>
#include <L1/Tool.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  //This needs some cleaning up, but it is automatically updated by SVN
  std::string Id = "$Id$";
  std::cout << "Executable tag: " << Id.substr(1,Id.length()-2) << std::endl;
  Core::Field< QCD::Gauge > myfield(8, 8);
  Tool::IO::load(&myfield, "../test/conf.88", Tool::IO::fileILDG);
  
  std::cout << myfield[0][0] << std::endl;
  return 0;
}
