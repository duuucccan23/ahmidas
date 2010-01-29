// $Id$
#include <string>
#include <iomanip>
#include <iostream>

#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Tool.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  //This needs some cleaning up, but it is automatically updated by SVN
  std::string Id = "$Id$";
  std::cout << "Executable tag: " << Id.substr(1,Id.length()-2) << std::endl;

  std::cout << "Reading gauge file conf.88 from test directory.\n";
  Core::Field< QCD::Gauge > myfield(8,8);
  Tool::IO::load(&myfield, "../../test/conf.88", Tool::IO::fileILDG);

  double plaqs = Tool::spatialPlaquette(myfield);
  double plaqt = Tool::temporalPlaquette(myfield);
  std::cout << "Summmed plaquette value: " << std::scientific << 0.5 * (plaqt + plaqs) << std::endl;
  std::cout << "\nThe summed plaquette value should match the one reported by viewing lime_contents on this configuration file.\n\n";
  
  std::cout << "Writing to conf.copy.88 in test directory.\n";
  Tool::IO::saveILDG(myfield, "../../test/conf.copy.88");
  std::cout << "To compare the binary info, extract the matching records from the two files with lime_extract_record, and diff.\n";
  return 0;
}
