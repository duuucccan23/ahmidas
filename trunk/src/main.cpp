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
  QCD::Tensor otherTens;
  for (size_t idx = 0; idx < 144; ++idx)
  {
    myTens[idx] = std::complex< double >(idx / 12, 0);
    otherTens[idx] = std::complex< double >(0, idx / 12);
  }
  myTens.rightMultiply(otherTens);
  for (size_t idx = 0; idx < 144; ++idx)
  {
    if (idx && idx % 12 == 0)
      std::cout << std::endl;
    std::cout << myTens[idx] << "\t";
  }
  std::cout << std::endl;
  return 0;
}
