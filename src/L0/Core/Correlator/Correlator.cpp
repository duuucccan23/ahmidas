#include "Correlator.ih"

template< >
Core::Field< int * > * Core::Correlator< Dirac::Matrix >::s_xRelative = NULL;

template< >
Core::Field< int * > * Core::Correlator< std::complex< double > >::s_xRelative = NULL;
