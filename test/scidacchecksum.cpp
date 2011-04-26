#include <iostream>
#include <L0/Ahmidas.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L1/Tool/IO.h>
#include <L0/Print.h>

int main(int argc, char **argv)
{
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <ildg_conf> <L> <T>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  Ahmidas my_ahmidas(&argc, &argv);
  Core::Field< QCD::Gauge > mygauge(atoi(argv[2]),atoi(argv[3]));
  Tool::IO::loadILDG(&mygauge, argv[1]);
  return 0;
}
