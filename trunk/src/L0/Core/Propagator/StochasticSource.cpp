#include "Propagator.ih"

namespace Core
{
  template<  >
  Propagator StochasticSource< 4 >::operator*(StochasticPropagator< 4 > const &sPropagator) const
  {

    Propagator tmp(sPropagator);
    assert (dynamic_cast< Propagator const *>(this)->T()==sPropagator.T() && dynamic_cast< Propagator const *>(this)->L()==sPropagator.L());

//     if (!(NComp == 4))
//     {
//       std::cerr << "Propagator StochasticSource< NComp >::operator*(StochasticPropagator< NComp > const &) const\n"
//           << "has not been implemented yet for NComp != 4" << std::endl;
//       exit(1);
//     }

    Propagator::iterator Itmp(tmp.begin()); //isolate() is called automatically here
    Propagator::const_iterator Is(begin());
    Propagator::const_iterator Ip(sPropagator.begin());

    while(Is != end())
    {
      std::transform(&((*Is)[0]), &((*Is)[0]) + 144, &((*Ip)[0]), &((*Itmp)[0]),
                     std::multiplies< std::complex< double > >());
      ++Itmp;
      ++Is;
      ++Ip;
    }
    assert(Itmp == tmp.end());
    assert(Ip == sPropagator.end());
    return tmp;
  }


  template< >
  Propagator StochasticSource< 12 >::createStochasticPropagator_fixedSink(StochasticPropagator< 12 > const &sPropagator, size_t const *sink) const
  {
    assert (dynamic_cast< Propagator const *>(this)->T()==sPropagator.T() &&
            dynamic_cast< Propagator const *>(this)->L()==sPropagator.L());
    Propagator tmp(dynamic_cast< Propagator const *>(this)->L(), dynamic_cast< Propagator const *>(this)->T());

    Propagator::iterator Itmp(tmp.begin()); //isolate() is called automatically here
    Propagator::const_iterator Is(begin());

    QCD::Tensor phi_sink(sPropagator(sink));

    while(Is != end())
    {
      std::transform(&((phi_sink)[0])      , &((phi_sink)[0]) +  12, &((*Itmp)[0]),
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[  0])));
      std::transform(&((phi_sink)[0]) +  12, &((phi_sink)[0]) +  24, &((*Itmp)[0]) +  12,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[ 13])));
      std::transform(&((phi_sink)[0]) +  24, &((phi_sink)[0]) +  36, &((*Itmp)[0]) +  24,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[ 26])));
      std::transform(&((phi_sink)[0]) +  36, &((phi_sink)[0]) +  48, &((*Itmp)[0]) +  36,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[ 39])));
      std::transform(&((phi_sink)[0]) +  48, &((phi_sink)[0]) +  60, &((*Itmp)[0]) +  48,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[ 52])));
      std::transform(&((phi_sink)[0]) +  60, &((phi_sink)[0]) +  72, &((*Itmp)[0]) +  60,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[ 65])));
      std::transform(&((phi_sink)[0]) +  72, &((phi_sink)[0]) +  84, &((*Itmp)[0]) +  72,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[ 78])));
      std::transform(&((phi_sink)[0]) +  84, &((phi_sink)[0]) +  96, &((*Itmp)[0]) +  84,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[ 91])));
      std::transform(&((phi_sink)[0]) +  96, &((phi_sink)[0]) + 108, &((*Itmp)[0]) +  96,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[104])));
      std::transform(&((phi_sink)[0]) + 108, &((phi_sink)[0]) + 120, &((*Itmp)[0]) + 108,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[117])));
      std::transform(&((phi_sink)[0]) + 120, &((phi_sink)[0]) + 132, &((*Itmp)[0]) + 120,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[130])));
      std::transform(&((phi_sink)[0]) + 132, &((phi_sink)[0]) + 144, &((*Itmp)[0]) + 132,
                     std::bind1st(std::multiplies< std::complex< double > >(), conj((*Is)[143])));
      ++Itmp;
      ++Is;
    }
    return tmp;
  }


  template< >
  Propagator StochasticSource< 1 >::createStochasticPropagator_fixedSink(StochasticPropagator< 1 > const &sPropagator, size_t const *sink) const
  {
    QCD::Tensor phi_sink(sPropagator(sink));
    assert (dynamic_cast< Propagator const *>(this)->T()==sPropagator.T() &&
            dynamic_cast< Propagator const *>(this)->L()==sPropagator.L());
    Propagator tmp(dynamic_cast< Propagator const *>(this)->L(), dynamic_cast< Propagator const *>(this)->T());

    Propagator::iterator Itmp(tmp.begin()); //isolate() is called automatically here
    Propagator::const_iterator Is(begin());

    while(Is != end())
    {
      std::transform(&((*Is)[0]), &((*Is)[0]) + 144, &(phi_sink[0]), &((*Itmp)[0]),
                     std::multiplies< std::complex< double > >());
      ++Itmp;
      ++Is;
    }
    return tmp;
  }

}
