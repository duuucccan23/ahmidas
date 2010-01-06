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

#pragma once

#include <algorithm>
#include <iostream>

#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>

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
    size_t                *d_references;

    Base::Weave< L, T >   &d_weave;
    Element               *d_field;
    size_t                *d_offsets;

    public:
      Field();
      Field(Element const &value);
      Field(Field< Element, L, T> const &other);

      template< typename Super >
      Field(Component< Super, L, T, Element> const &component); //Creation of a field through a Component of another field.

      template< typename Super >
      Field(hcComponent< Super, L, T, Element> const &component); //Creation of a field through a hcComponent of another field.

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
      Component< Element, L, T, Atom > component(size_t const component);

      Field< Element, L, T > &shift(Base::SpaceTimeIndex idx, Base::Direction shift);

#include "Field/Field.operators"

      hcField< Element, L, T > dagger() const;

      // Unsafe memory access: ONLY when you really control all repercussions
      Element * const raw();

      // Indexing
      Element &operator[](size_t const idx);
      Element const &operator[](size_t const idx) const;

      Element &physicalIndex(size_t const idx);
      Element &fastPhysicalIndex(size_t const idx);
      Element const &physicalIndex(size_t const idx) const;
      Element const &constPhysicalIndex(size_t const idx) const;

      Element &memoryIndex(size_t const idx);
      Element &fastMemoryIndex(size_t const idx);
      Element const &memoryIndex(size_t const idx) const;
      Element const &constMemoryIndex(size_t const idx) const;

      Element *at(size_t x, size_t y, size_t z, size_t t);
      Element *fastAt(size_t x, size_t y, size_t z, size_t t);
      Element const *at(size_t x, size_t y, size_t z, size_t t) const;
      Element const *constAt(size_t x, size_t y, size_t z, size_t t) const;

      size_t size() const; //Enumerates the number of Elements in this field (for looping)
      size_t spatialSize() const; // (for looping over a single timeslice)
      void isolate();
      void fill(Element const &element); // Flush a field with a constant quantity

    private:
      size_t shiftIdxToZero(size_t const idx) const;
      size_t shiftIdxToOffset(size_t const idx) const;

      void destroy();
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

      size_t size() const;
  };
}

#include "Field/Field.inlines"
#include "Field/Field.iterator.inlines"
#include "Field/Field.const_iterator.inlines"
#include "Field/Field_Field_a.template"
#include "Field/Field_Field_b.template"
#include "Field/Field_Field_c.template"
#include "Field/Field_Field_d.template"
#include "Field/Field_Field_e.template"
#include "Field/Field_Field_f.template"
#include "Field/Field_destroy.template"
#include "Field/Field_operator_eq.template"
#include "Field/Field_shift.template"
#include "Field/Field_isolate.template"
#include "Field/Field.operators.templates"

#include "Field/hcField.inlines"
