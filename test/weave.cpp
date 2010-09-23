#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>
#include <L0/Debug.h>

int main(int argc, char **argv)
{
  std::ostringstream ostr("Weaving\n", std::ios::ate);
  Base::Weave myWeave = Base::Weave(4, 8);
  ostr << "LocalVolume() " << myWeave.localVolume() << std::endl << myWeave.dim(Base::idx_X) << '\t'
       << myWeave.dim(Base::idx_Y) << '\t' << myWeave.dim(Base::idx_Z) << '\t'
       << myWeave.dim(Base::idx_T) << std::endl << myWeave.localSize(Base::idx_X) << '\t'
       << myWeave.localSize(Base::idx_Y) << '\t' << myWeave.localSize(Base::idx_Z) << '\t'
       << myWeave.localSize(Base::idx_T);
  Debug(ostr.str(), std::cout);
  return 0;
}
