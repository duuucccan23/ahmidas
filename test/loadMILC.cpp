#include <iostream>
#include <L0/Ahmidas.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Tool.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);
  Core::Field< QCD::Gauge > field(16,6);
  Tool::IO::load(&field, "../../test/lat.sample.16x6.milc", Tool::IO::fileMILC);
  double sp = Tool::spatialPlaquette(field);
  double tp = Tool::temporalPlaquette(field);
  std::cout << "Spatial plaquette:  " << sp << ", should be 0.396411.\n";
  std::cout << "Temporal plaquette: " << tp << ", should be 0.399299.\n";

  if ((fabs(sp - 0.396411) < 1e-6) && (fabs(tp - 0.399299) < 1e-6))
    return 0;
  return 1;
}
