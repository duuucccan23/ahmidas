/*******************************************************************************
 *    This file is part of Ahmidas.                                            *
 *                                                                             *
 *    Ahmidas is free software: you can redistribute it and/or modify          *
 *    it under the terms of the GNU General Public License as published by     *
 *    the Free Software Foundation, either version 3 of the License, or        *
 *    (at your option) any later version.                                      *
 *                                                                             *
 *    Ahmidas is distributed in the hope that it will be useful,               *
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *    GNU General Public License for more details. You should have received    *
 *    a copy of the GNU General Public License along with Ahmidas. If not,     *
 *    see <http://www.gnu.org/licenses/>.                                      *
 *                                                                             *
 *    Although we hereby grant the right to use Ahmidas or its parts as you    *
 *    see fit, we kindly ask that you mention its use in any scientific        *
 *    papers that benefit from it. In that case, please refer to               *
 *    < add reference here >.                                                  *
 *******************************************************************************/

#ifndef GUARD_CORE_FIELD_H
#define GUARD_CORE_FIELD_H

#include <algorithm>
#include <L0/Core/Unifold.h>

/* A field is a 4 dimensional object with an X, Y, Z and T dimension. X = Y = Z = L.
 * 
 */

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
    size_t           *d_references;

    Unifold< L, T >  &d_unifold;
    Element          *d_field;
    size_t           *d_offsets;

    public:
      Field();
      Field(Element const &value);
      Field(Field const &other);

      explicit Field(hcField< Element, L, T > const &other);
      Field< Element, L, T > &operator=(Field< Element, L, T > const &other);

      ~Field();

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

      Element &getMemoryIndex(size_t const idx);
      Element const &getMemoryIndex(size_t const idx) const;

    private:
      size_t shiftIdxToZero(size_t const idx) const;
      size_t shiftIdxToOffset(size_t const idx) const;

      void destroy();
      void isolate();
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
#include "Field/Field_Field_b.template"
#include "Field/Field_Field_c.template"
#include "Field/Field_Field_d.template"
#include "Field/Field_~Field.template"
#include "Field/Field_operator_eq.template"
#include "Field/Field_shift.template"
#include "Field/Field_isolate.template"
#include "Field/Field_destroy.template"

#include "Field/hcField.inlines"

#endif