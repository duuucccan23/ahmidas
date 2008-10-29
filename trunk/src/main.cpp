#include <iomanip>
#include <iostream>

#include <L0/Base/Base.h>
#include <L0/Base/IO.h>
#include <L0/Base/Random.h>
#include <L0/Core/Field.h>
#include <L0/Core/TMatrix.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Tensor.h>
#include <L1/Source/Point.h>
#include <L1/Tool.h>
#include <L2/Contract/Wick.h>
#include <L2/Contract/TwoPoint.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Spinor, 24, 48 > source;
//  std::cout << "Raw mem source:\n" << source[0][0] << source[0][1] << source[0][2] << source[0][3] << "\n" << source[1][0] << "\n";
  Tool::randomize(&source);
//  std::cout << "Zeroed source:\n" << source[0][0] << source[0][1] << source[0][2] << source[0][3] << "\n" << source[1][0] << "\n";
  Base::IO::saveScidac(source, "../test/ahmsource0.2448");
  Base::IO::saveScidac(source, "../test/ahmsource1.2448");
  Base::IO::saveScidac(source, "../test/ahmsource2.2448");
  Base::IO::saveScidac(source, "../test/ahmsource3.2448");
  /*
  Core::Field< QCD::Spinor, 24, 48 > source;
  Base::IO::loadScidac(&source, "../test/smearsource0.2448");
  std::cout << "First source:\n" << source[0][0] << source[0][1] << source[0][2] << source[0][3] << "\n" << source[1][0] << "\n";
  Base::IO::loadScidac(&source, "../test/smearsource1.2448");
  std::cout << "Second source:\n" << source[0][0] << source[0][1] << source[0][2] << source[0][3] << "\n" << source[1][0] << "\n";
  Base::IO::loadScidac(&source, "../test/smearsource2.2448");
  std::cout << "Third source:\n" << source[0][0] << source[0][1] << source[0][2] << source[0][3] << "\n" << source[1][0] << "\n";
  Base::IO::loadScidac(&source, "../test/smearsource3.2448");
  std::cout << "Last source:\n" << source[0][0] << source[0][1] << source[0][2] << source[0][3] << "\n" << source[1][0] << "\n";
  */
  /*Core::Propagator< 24, 48 > propagator;
  for (Core::Propagator< 24, 48 >::iterator_full iter = propagator.begin(); iter != propagator.end(); ++iter)
    Base::IO::loadScidac(&(*iter), "../test/inv.2448");

  Source::Point< 24, 48 > pointSource(3, 2, 1);

  Core::TMatrix< 24, 48 > tmatrix = Contract::wick(pointSource, propagator);
  Core::Correlator< 24, 48 > correlator = Contract::twoPoint(tmatrix.dagger(), tmatrix);

  for (size_t idx = 0; idx < 48; ++idx)
    std::cout << std::setw(2) << idx << "\t" << correlator[idx].real() << std::endl;
    */  
  return 0;
}
