#include "Propagator.ih"

//multply only on sink by (1+-g5)/sqrt(2)

namespace Core
{
  template< >
  void StochasticPropagator <1> ::rotateToPhysicalBasis(bool const sign)
  {
    isolate();

    Dirac::Gamma< 5 > gamma5;

    StochasticPropagator<1> tmp(*this);
    tmp.rightMultiply(gamma5);
    if(sign==0) tmp *= std::complex< double >(0,-1);
    else tmp *= std::complex< double >(0,+1);
    (*this) += tmp;

    (*this) *= 1/sqrt(2.0);
  }

  template< >
  void StochasticPropagator <4>::rotateToPhysicalBasis(bool const sign)
  {
    isolate();

    Dirac::Gamma< 5 > gamma5;

    StochasticPropagator<4> tmp(*this);
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

/*
  template <>
  void StochasticPropagator<12>::rotateToPhysicalBasis(bool const sign)
  {
    isolate();

    Dirac::Gamma< 5 > gamma5;

    StochasticPropagator<12> tmp(*this);
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
*/
}
