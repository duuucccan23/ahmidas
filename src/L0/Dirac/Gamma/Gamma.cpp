#include "Gamma.ih"

namespace Dirac
{
  namespace
  {
    std::complex< double > const I = std::complex< double >(0, 1);
  }


  // identity
  template< >
  size_t const Gamma< -1 >::s_perm[4] = {0, 1, 2, 3};
  template< >
  std::complex< double > const Gamma< -1 >::s_sign[4] = {1, 1, 1, 1};

  // Single Gamma terms
  template< >
  size_t const Gamma< 1 >::s_perm[4] = {3, 2, 1, 0};
  template< >
  std::complex< double > const Gamma< 1 >::s_sign[4] = {-I, -I, I, I};

  template< >
  size_t const Gamma< 2 >::s_perm[4] = {3, 2, 1, 0};
  template< >
  std::complex< double > const Gamma< 2 >::s_sign[4] = {1, -1, -1, 1};

  template< >
  size_t const Gamma< 3 >::s_perm[4] = {2, 3, 0, 1};
  template< >
  std::complex< double > const Gamma< 3 >::s_sign[4] = {-I, I, I, -I};

  template< >
  size_t const Gamma< 4 >::s_perm[4] = {2, 3, 0, 1};
  template< >
  std::complex< double > const Gamma< 4 >::s_sign[4] = {-1, -1, -1, -1};

  template< >
  size_t const Gamma< 5 >::s_perm[4] = {0, 1, 2, 3};
  template< >
  std::complex< double > const Gamma< 5 >::s_sign[4] = {1, 1, -1, -1};

  // Pseudoscalar combinations
  template< >
  size_t const Gamma< 15 >::s_perm[4] = {3, 2, 1, 0};
  template< >
  std::complex< double > const Gamma< 15 >::s_sign[4] = {I, I, I, I};

  template< >
  size_t const Gamma< 25 >::s_perm[4] = {3, 2, 1, 0};
  template< >
  std::complex< double > const Gamma< 25 >::s_sign[4] = {1, -1, 1, -1};

  template< >
  size_t const Gamma< 35 >::s_perm[4] = {2, 3, 0, 1};
  template< >
  std::complex< double > const Gamma< 35 >::s_sign[4] = {I, -I, I, -I};

  template< >
  size_t const Gamma< 45 >::s_perm[4] = {2, 3, 0, 1};
  template< >
  std::complex< double > const Gamma< 45 >::s_sign[4] = {1, 1, -1, -1};

  template< >
  size_t const Gamma< 51 >::s_perm[4] = {3, 2, 1, 0};
  template< >
  std::complex< double > const Gamma< 51 >::s_sign[4] = {-I, -I, -I, -I};

  template< >
  size_t const Gamma< 52 >::s_perm[4] = {3, 2, 1, 0};
  template< >
  std::complex< double > const Gamma< 52 >::s_sign[4] = {-1, 1, -1, 1};

  template< >
  size_t const Gamma< 53 >::s_perm[4] = {2, 3, 0, 1};
  template< >
  std::complex< double > const Gamma< 53 >::s_sign[4] = {-I, I, -I, I};

  template< >
  size_t const Gamma< 54 >::s_perm[4] = {2, 3, 0, 1};
  template< >
  std::complex< double > const Gamma< 54 >::s_sign[4] = {-1, -1, 1, 1};
}
