#include <iostream>

#include <L0/Base/Base.h>
#include <L0/Base/IO.h>
#include <L0/Base/Random.h>
#include <L0/Core/Field.h>
#include <L0/Core/TMatrix.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Tensor.h>
#include <L1/Source/Stochastic.h>
#include <L1/Source/Full.h>
#include <L1/Tool.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Spinor, 24, 48 > testspinor = Base::IO::loadScidac< QCD::Spinor, 24, 48 >("../test/smearsource0.2448");
  std::cout << testspinor[0][0];
  Base::IO::saveILDG(testspinor, "../test/testspinor.2448");
  Source::Full< 8, 8 > mysource = Source::Full< 8, 8 >();
  Tool::IO::saveScidac(mysource, "../test/fulltest.88");
  return 0;
}
