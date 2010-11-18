#include <iostream>
#include <L0/Base/Weave.h>
#include <L1/Tool/IO/Lime/Reader.h>

int main(int argc, char **argv)
{
  Base::Weave weave(4, 4);
  Tool::IO::Lime::Reader limeTest(&weave, argv[1]);
  // Something more intelligent should be done here, probably.
  // Maybe we should construct a special test file for LIME.
  // This particular file is in fact ILDG and needs more care.
  return 0;
}
