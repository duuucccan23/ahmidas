#include <iostream>
#include <L0/QCD/Tensor.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>
#include <L0/Base/Base.h>
#include <L0/Base/IO.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Spinor, 24, 48 > myspinor;
  Base::IO::loadScidac(&myspinor, "../test/inv.2448");

  QCD::Tensor myTens;
  for (size_t idx = 0; idx < 144; ++idx)
    myTens[idx] = std::complex< double >(idx % 12, idx / 12);
  myTens.leftMultiply(myTens.dagger());
  for (size_t idx = 0; idx < 144; ++idx)
  {
    if (idx != 0 && ((idx % 12) == 0))
      std::cout << std::endl;
    std::cout << myTens[idx] << " ";
  }
  std::cout << std::endl;
  return 0;
}
