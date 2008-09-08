#ifndef GUARD_CORE_FIELD_H
#define GUARD_CORE_FIELD_H

#include <cassert>
#include <complex>
#include <L0/IO/Lime/Reader.h>

#include <mpi.h>
#include <L0/Core/Grid.h>
#include <L0/Core/Com.h>

/* A field is a 4 dimensional object, with an X,Y,Z and T dimension. X=Y=Z=L.
 *
 */

namespace IO
{
  template< typename IOClass, typename Element, size_t L, size_t T >
  class Bridge;
}

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
    template< typename fIOClass, typename fElement, size_t fL, size_t fT >
    friend class IO::Bridge;

    size_t               *d_references;

    Com< Element, L, T>  &d_com;
    Element              *d_field;
    size_t               *d_offsets;

    public:
      Field();
      Field(Element const &value);
      Field(Field const &other);

      explicit Field(hcField< Element, L, T > const &other);
      Field< Element, L, T > &operator=(Field< Element, L, T > const &other);

      ~Field();

      template< typename Precision >
      void readFromFile(char const* fileName, char const* fileType);

      void increaseIdx(size_t *idx) const;
      void decreaseIdx(size_t *idx) const;

      Element &element(size_t const *idx);
      Element const &element(size_t const *idx) const;

#include "Field/Field.iterator"
#include "Field/Field.const_iterator"

      iterator begin();
      iterator end();

      const_iterator begin() const;
      const_iterator end() const;

      template< typename Atom >
      Component< Element, L, T, Atom > component(SpaceTimeIndex idx);

      Field< Element, L, T > &shift(SpaceTimeIndex idx, Direction shift);

#include "Field/Field.operators"

      hcField< Element, L, T > dagger() const;

      void averageTimeSlice(std::complex< double > *result);

      template< typename IOClass >
      void loadDataFromIO(IOClass &inputIO);

    private:
      size_t moveBufferToData(char *fileBuffer, size_t written, size_t precision);
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

}

#include "Field/Field.inlines"
#include "Field/Field.iterator.inlines"
#include "Field/Field.const_iterator.inlines"
#include "Field/Field_Field_a.template"
#include "Field/Field_Field_c.template"
#include "Field/Field_Field_d.template"
#include "Field/Field_Field_e.template"
#include "Field/Field_~Field.template"
#include "Field/Field_decreaseIndex.template"
#include "Field/Field_increaseIndex.template"
#include "Field/Field_loadDataFromIO.template"
#include "Field/Field_moveBufferToData.template"
#include "Field/Field_operator_eq.template"
#include "Field/Field_shift.template"

#endif
