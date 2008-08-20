#ifndef GUARD_CONTRACT_EXPRESSION_H
#define GUARD_CONTRACT_EXPRESSION_H

namespace Contract
{
  template< typename FieldLeft, typename FieldRight >
  class Expression
  {
    FieldLeft &d_left;
    FieldRight &d_right;

    public:
      Expression(FieldLeft &left, FieldRight &right);
  };

  template< size_t Index, size_t L, size_t T >
  Expression< Dirac::Gamma< Index >, Core::Field< QCD::Spinor, L, T > >
    operator*(Dirac::Gamma< Index > gamma, Core::Field< QCD::Spinor, L, T > field);
}

#include "Expression/Expression.inlines"

#endif
