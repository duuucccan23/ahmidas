#include <iostream>
#include <L0/Core/Correlator.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Source/Point.h>
#include <L2/Contract/Wick.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Spinor, 8, 8 > myfield;
  for (size_t idx = 0; idx < 64 * 64; ++idx)
    myfield[idx].setToRandom();
  Source::Point< 8, 8 > mysource(1, 1, 1);
  Core::Correlator< 8, 8 > mycorr = Contract::wick(mysource, Base::gam_1, Base::col_RED, myfield);
  return 0;
}
