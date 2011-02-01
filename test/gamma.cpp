
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

  Gamma< 1 > gamma1;
  Gamma< 2 > gamma2;
  Gamma< 3 > gamma3;
  Gamma< 4 > gamma4;
  Gamma< 5 > gamma5;
  complex <double >  plus_i(0,1);
  complex <double >  minus_i(0,-1);


  QCD::Spinor spinor1(refData);
  spinor1.leftMultiply(gamma1);

  for (int idx=0; idx<3; idx++)
  {
	assert(spinor1[0][idx] == minus_i*refData[3*3+idx]);
    assert(spinor1[1][idx] == minus_i*refData[2*3+idx]);
    assert(spinor1[2][idx] == plus_i*refData[1*3+idx]);
    assert(spinor1[3][idx] == plus_i*refData[0*3+idx]);
  }

  Print("left multiplication by gamma1 ok");


  QCD::Spinor spinor2(refData);
  spinor2.leftMultiply(gamma2);

  for (int idx=0; idx<3; idx++)
  {
	assert(spinor2[0][idx] == -refData[3*3+idx]);
    assert(spinor2[1][idx] ==  refData[2*3+idx]);
    assert(spinor2[2][idx] ==  refData[1*3+idx]);
    assert(spinor2[3][idx] == -refData[0*3+idx]);
  }

  Print("left multiplication by gamma2 ok");


  QCD::Spinor spinor3(refData);
  spinor3.leftMultiply(gamma3);

  for (int idx=0; idx<3; idx++)
  {
	assert(spinor3[0][idx] == minus_i*refData[2*3+idx]);
    assert(spinor3[1][idx] == plus_i*refData[3*3+idx]);
    assert(spinor3[2][idx] == plus_i*refData[0*3+idx]);
    assert(spinor3[3][idx] == minus_i*refData[1*3+idx]);
  }

  Print("left multiplication by gamma3 ok");


  QCD::Spinor spinor4(refData);
  spinor4.leftMultiply(gamma4);

  for (int idx=0; idx<3; idx++)
  {
	assert(spinor4[0][idx] == -refData[2*3+idx]);
    assert(spinor4[1][idx] == -refData[3*3+idx]);
    assert(spinor4[2][idx] == -refData[0*3+idx]);
    assert(spinor4[3][idx] == -refData[1*3+idx]);
  }

  Print("left multiplication by gamma4 ok");



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

  Print("left multiplication by gamma5 ok");

  QCD::Spinor spinor6(refData);
  spinor6.leftMultiply(gamma5);
  spinor6.leftMultiply(gamma1);


  QCD::Spinor spinor7(refData);
  spinor7.leftMultiply(gamma1);
  spinor7.leftMultiply(gamma5);



  // loop over color entries
  for (int idx=0; idx<3; idx++)
  {
    assert(spinor6[0][idx] + spinor7[0][idx] ==  0.);
    assert(spinor6[1][idx] + spinor7[1][idx] ==  0.);
    assert(spinor6[2][idx] + spinor7[2][idx] ==  0.);
    assert(spinor6[3][idx] + spinor7[3][idx] ==  0.);
  }

  Print("left multiplication by {gamma5 ,gamma1} ok");






  return EXIT_SUCCESS;
}

