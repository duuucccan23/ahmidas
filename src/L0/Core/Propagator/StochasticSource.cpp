#include "Propagator.ih"

namespace Core
{
  template< size_t NComp >
  Propagator StochasticSource< NComp >::operator*(StochasticPropagator< NComp > const &sPropagator) const
  {

    std::cerr << "Propagator StochasticSource< NComp >::operator*(StochasticPropagator< NComp > const &) const\n"
              << "has not been implemented yet" << std::endl;
    exit(1);

    Propagator tmp(sPropagator);
//     //assert (T()==propagator.T() && L()==propagator.L());
// 
// 
// 
//     Propagator::iterator Itmp(tmp.begin()); //isolate() is called automatically here
//     Propagator::const_iterator Is(begin());
//     Propagator::const_iterator Ip(begin());
// 
//     while(Is != end())
//     {
//       std::transform(&((*Is)[0]), &((*Is)[0]) + 144, &((*Ip)[0]), &((*Itmp)[0]),
//                      std::multiplies< std::complex< double > >());
//       ++Itmp;
//       ++Is;
//       ++Ip;
//     }
//     assert(Itmp == tmp.end());
//     assert(Ip == propagator.end());
    return tmp;
  }

}
