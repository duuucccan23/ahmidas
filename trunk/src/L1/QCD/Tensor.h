#ifndef GUARD_QCD_TENSOR_H
#define GUARD_QCD_TENSOR_H

#include <L0/Core/Core.h>
#include <L0/Core/Field.h>
#include <L0/Core/Grid.h>
#include <L0/QCD/Spinor.h>

namespace QCD
{
  template< size_t L, size_t T >
  class Tensor
  {
    Core::Field< Spinor, L, T > *d_components[12];

    public:
      Tensor(Core::Grid< L, T > &grid);
      Tensor(Tensor const &other);
      Tensor(Core::Field< Spinor, L, T > *data);
      ~Tensor();

#include "Tensor/Tensor.iterator"
#include "Tensor/Tensor.const_iterator"

      iterator_full begin();
      iterator_full end();

      const_iterator_full begin() const;
      const_iterator_full end() const;

      iterator_colour begin(Core::ColourIndex const idx);
      iterator_colour end(Core::ColourIndex const idx);

      const_iterator_colour begin(Core::ColourIndex const idx) const;
      const_iterator_colour end(Core::ColourIndex const idx) const;

      iterator_dirac begin(Core::DiracIndex const idx);
      iterator_dirac end(Core::DiracIndex const idx);

      const_iterator_dirac begin(Core::DiracIndex const idx) const;
      const_iterator_dirac end(Core::DiracIndex const idx) const;
  };
}

#include "Tensor/Tensor.inlines"
#include "Tensor/Tensor.iterator.inlines"
#include "Tensor/Tensor.const_iterator.inlines"

#endif
