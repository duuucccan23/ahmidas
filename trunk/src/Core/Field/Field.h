#ifndef GUARD_CORE_FIELD_H
#define GUARD_CORE_FIELD_H

#include <cassert>
#include <complex>
#include <mpi.h>
#include <iostream>
#include <Core/Grid/Grid.h>
#include <Lime/Reader/Reader.h>

namespace Core
{
  template< typename Element, size_t L, size_t T, typename Atom >
  class Component;

  template< typename Element, size_t L, size_t T, typename Atom >
  class hcComponent;


  template< typename Element, size_t L, size_t T >
  class hcField;

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

      explicit Field(hcField< Element, L, T > const &other);
      Field< Element, L, T > &operator=(Field< Element, L, T > const &other);

      ~Field();

      template< typename Precision >
      void readFromFile(char const* fileName, char const* fileType);

      Core::Grid< L, T > const &grid() const;

      void increaseIdx(size_t *idx) const;
      void decreaseIdx(size_t *idx) const;

      Element &element(size_t const *idx);

#include "Field.iterator"
#include "Field.const_iterator"

      iterator begin();
      iterator end();

      const_iterator begin() const;
      const_iterator end() const;

      template< typename Atom >
      Component< Element, L, T, Atom > component(SpaceTimeIndex idx);

      Field< Element, L, T > &shift(SpaceTimeIndex idx, Direction shift);

#include "Field.operators"

      hcField< Element, L, T > dagger() const;
      void reunitarize();

      void averageTimeSlice(std::complex< double > *result);

    private:
      void setSurfaces();
      size_t moveBufferToData(size_t written);
  };

  template< typename Element, size_t L, size_t T >
  class hcField
  {
    friend class Field< Element, L, T >;
    Field< Element, L, T > const &d_parent;

    hcField(Field< Element, L, T > const &parent);

    public:
      Field< Element, L, T > const &dagger() const;

      typename Field< Element, L, T >::const_iterator begin() const;
      typename Field< Element, L, T >::const_iterator end() const;
  };

  std::complex< double > plus(std::complex< double > const &left, std::complex< double > const &right);
}

#include "Field.inlines"
#include "Field.iterator.inlines"
#include "Field.const_iterator.inlines"
#include "Field_a.template"
#include "Field_b.template"
#include "Field_c.template"
#include "Field_d.template"
#include "Field_e.template"
#include "Field_f.template"
#include "Field_g.template"
#include "Field_h.template"
#include "Field_i.template"
#include "Field_j.template"
#include "Field_k.template"

#endif
