#include "Propagator.ih"

//multply by left and right by (1+-g5)/sqrt(2)

namespace Core
{
  void Propagator::rotateToPhysicalBasis(bool const sign)
  {
    isolate();

    Dirac::Gamma< 5 > gamma5;

    Propagator tmp(*this);
    tmp *= gamma5;
    if(sign==0) tmp *= std::complex< double >(0,-1);
    else tmp *= std::complex< double >(0,+1);
    (*this) += tmp;

    tmp = (*this);
    tmp.rightMultiply(gamma5);
    if(sign==0) tmp *= std::complex< double >(0,-1);
    else tmp *= std::complex< double >(0,+1);
    (*this) += tmp;

    (*this) *= 0.5; //this is (1/sqrt(2))^2
  }
}
