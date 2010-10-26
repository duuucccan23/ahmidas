#include "Propagator.ih"

//multply by left and right by (1+-g5)/sqrt(2)

namespace Core
{
  void Propagator::rotateToPhysicalBasis(bool const sign)
  {
    isolate();

    std::complex< double > const Ifactor = sign ?  std::complex< double >(0,-1) : std::complex< double >(0,+1);
    Dirac::Gamma< 5 > gamma5;

    Propagator tmp(*this);
    tmp *= gamma5;
    tmp *= Ifactor;
    (*this) += tmp;

    tmp = (*this);
    tmp.rightMultiply(gamma5);
    tmp *= Ifactor;
    (*this) += tmp;

    (*this) *= 0.5; //this is (1/sqrt(2))^2
  }
}
