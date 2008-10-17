#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Source/Point.h>
#include <L2/Contract/Wick.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.88");
  Source::Point< 8, 8 > mysource(1, 1, 1);
  Core::Correlator< 8 > mycorr = Contract::wick(mysource, Base::col_RED, Base::dir_0, myfield);
  return 0;
}
