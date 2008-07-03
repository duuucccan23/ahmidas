#ifndef GUARD_FIELDS_BUFFER_H
#define GUARD_FIELDS_BUFFER_H

#include "../../Core/Grid/Grid.h"

namespace Core
{
  template< typename Element >
  class Buffer
  {
    size_t   d_size;
    Element *d_data;

    public:
      template< size_t L, size_t T >
      Buffer(Core::Grid< L, T > const &grid);

      template< size_t L, size_t T >
      Buffer(Core::Grid< L, T > const &grid, Element const &defValue);

      Buffer(Buffer< Element > &other); // Will pilfer resources!
      Buffer(Buffer< Element > &other, Element const &defValue); // Will pilfer resources!

      ~Buffer();

      template< typename OtherElement, size_t L, size_t T >
      void leftMultiply(Field< OtherElement, L, T > const &other);

      template< typename OtherElement, size_t L, size_t T >
      void leftMultiply(hcField< OtherElement, L, T > const &other);

      template< typename OtherElement, size_t L, size_t T, typename Atom >
      void leftMultiply(Component< OtherElement, L, T, Atom > const &other);

      template< typename OtherElement, size_t L, size_t T, typename Atom >
      void leftMultiply(hcComponent< OtherElement, L, T, Atom > const &other);

      template< typename OtherElement, size_t L, size_t T >
      void rightMultiply(Field< OtherElement, L, T > const &other);

      template< typename OtherElement, size_t L, size_t T >
      void rightMultiply(hcField< OtherElement, L, T > const &other);

      template< typename OtherElement, size_t L, size_t T, typename Atom >
      void rightMultiply(Component< OtherElement, L, T, Atom > const &other);

      template< typename OtherElement, size_t L, size_t T, typename Atom >
      void rightMultiply(hcComponent< OtherElement, L, T, Atom > const &other);

      typedef Element* iterator;
      iterator begin();
      iterator end();
  };

  // Functor definitions for a clean interface
  template< typename Element >
  class add_buffer
  {
    Element *d_pBuffer;

    public:
      add_buffer(Buffer< Element > &buffer);

      template< typename fType >
      void operator()(fType &field);
  };

  template< typename Element  >
  add_buffer< Element > add(Buffer< Element > &buffer); // For further convenience
}

#include "Buffer.inlines"
#include "Buffer.functor.inlines"

#endif
