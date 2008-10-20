#include <iostream>
#include <L0/QCD/Tensor.h>

int main(int argc, char **argv)
{
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
