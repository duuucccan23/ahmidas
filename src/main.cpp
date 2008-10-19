#include <iostream>
#include <L0/Core/Propagator.h>
#include <L0/Core/TMatrix.h>
#include <L1/Source/Point.h>
#include <L1/Tool.h>
#include <L2/Contract/Wick.h>

int main(int argc, char **argv)
{
  Core::Propagator< 12, 8 > myprop;
  Tool::randomize(&myprop);

  Source::Point< 12, 8 > mysink(3, 5, 3);

  Core::TMatrix< 12, 8 > mytmat = Contract::wick(mysink, myprop);

  return 0;
}
