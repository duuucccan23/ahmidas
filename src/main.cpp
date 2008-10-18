#include <iostream>
#include <L0/Core/Correlator.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Source/Point.h>
#include <L1/Tool.h>
#include <L2/Contract/Wick.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Spinor, 8, 16 > myfield;
  Source::Point< 8, 16 > mysource(1, 1, 1);
  Tool::randomize(&myfield);
  Core::Correlator< 8, 16 > mycorr = Contract::wick(mysource, Base::gam_1, Base::col_RED, myfield);
  std::cout << "Result of correlator contraction of point sink on random spinor field:" << std::endl;
  for (size_t idx = 0; idx < 16; ++idx)
    std::cout << "T = " << idx << "\t" << mycorr[idx] << std::endl;
  return 0;
}
