#ifndef GUARD_CORE_FIELD_H
#define GUARD_CORE_FIELD_H

#include "../Grid/Grid.h"

namespace Core
{
  template< typename Element, typename Atom >
  class Component;

  template< typename Element, size_t L, size_t T >
  class Field
  {
    Grid< L, T >  &d_grid;
    size_t         d_offsets[4];
    size_t         d_bufferSize;
    MPI::Datatype  d_surfaces[4];

    Element  *d_buffer;
    Element  *d_field;

    public:
      Field(Core::Grid< L, T > &grid);
      Field(Field const &other);
      Field(Core::Grid< L, T > &grid, Element const &value);
      Field(Field const &other, Element const &value);
      ~Field();

      void readFromFile(char const* fileName); /// NOT YET IMPLEMENTED!
      Core::Grid< L, T > const &grid() const;

      void increaseIdx(short *idx) const;
      void decreaseIdx(short *idx) const;

      Element &element(short const *idx);

#include "Field.iterator"

      iterator begin();
      iterator end();

      template< typename Atom >
      Component< Element, Atom > component(size_t idx);

      Field< Element, L, T > &shift(SpaceTimeIndex idx, Direction shift);

#include "Field.operators"

      void reunitarize();
         
      void averageTimeSlice(std::complex< double > *result);

    private:
      void setSurfaces();
  };
}

#include "Field.inlines"
#include "Field.iterator.inlines"
#include "Field.templates"

#endif
