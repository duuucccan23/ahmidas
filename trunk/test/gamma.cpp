
#include <cassert>
#include <complex>

#include <L0/Print.h>
#include <L0/Ahmidas.h>
#include <L0/Dirac/Gamma.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Tensor.h>


using std::complex;
using Dirac::Gamma;

int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  Print("Testing some gamma multiplications...");

  complex< double > const refData[12] = {
      complex< double >( 1,  2), std::complex< double >( 3,  4),
      complex< double >( 5,  6), std::complex< double >( 7,  8),
      complex< double >( 9, 10), std::complex< double >(11, 12),
      complex< double >(13, 14), std::complex< double >(15, 16),
      complex< double >(17, 18), std::complex< double >(19, 20),
      complex< double >(21, 22), std::complex< double >(23, 24)};

  Gamma< 5 > gamma5;

  QCD::Spinor spinor(refData);
  spinor.leftMultiply(gamma5);

  // loop over color entries
  for (int idx=0; idx<3; idx++)
  {
    assert(spinor[0][idx] ==  refData[0*3+idx]);
    assert(spinor[1][idx] ==  refData[1*3+idx]);
    assert(spinor[2][idx] == -refData[2*3+idx]);
    assert(spinor[3][idx] == -refData[3*3+idx]);
  }

  Print("That seems to work.");
  return EXIT_SUCCESS;
}

