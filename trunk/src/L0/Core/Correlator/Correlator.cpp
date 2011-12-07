#include "Correlator.ih"

class complex12;
class complex16;
class complex96;
class complex192;

template<  >
Core::Field< int * > * Core::Correlator< complex192 >::s_xRelative = NULL;

template<  >
Core::Field< int * > * Core::Correlator< complex96 >::s_xRelative = NULL;

template<  >
Core::Field< int * > * Core::Correlator< complex16 >::s_xRelative = NULL;

template<  >
Core::Field< int * > * Core::Correlator< complex12 >::s_xRelative = NULL;

template< >
Core::Field< int * > * Core::Correlator< Dirac::Matrix >::s_xRelative = NULL;

template< >
Core::Field< int * > * Core::Correlator< std::complex< double > >::s_xRelative = NULL;

template<  >
Core::Field< size_t > * Core::Correlator< complex192 >::s_timelabel = NULL;

template<  >
Core::Field< size_t > * Core::Correlator< complex16 >::s_timelabel = NULL;

template<  >
Core::Field< size_t > * Core::Correlator< complex96 >::s_timelabel = NULL;

template<  >
Core::Field< size_t > * Core::Correlator< complex12 >::s_timelabel = NULL;

template< >
Core::Field< size_t > * Core::Correlator< Dirac::Matrix >::s_timelabel = NULL;

template< >
Core::Field< size_t > * Core::Correlator< std::complex< double > >::s_timelabel = NULL;

