#include "Correlator.ih"

class complex96;
class complex192;

template<  >
Core::Field< int * > * Core::Correlator< complex192 >::s_xRelative = NULL;

template<  >
Core::Field< int * > * Core::Correlator< complex96 >::s_xRelative = NULL;

template< >
Core::Field< int * > * Core::Correlator< Dirac::Matrix >::s_xRelative = NULL;

template< >
Core::Field< int * > * Core::Correlator< std::complex< double > >::s_xRelative = NULL;
