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
#include <complex>

#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>

namespace Core
{
  template< typename Element, typename Atom >
  class Component;

  template< typename Element, typename Atom >
  class hcComponent;

  template< typename Element >
  class hcField;

  template< typename Element >
  class Field
  {
    size_t                *d_references;

    Base::Weave            d_weave;
    Element               *d_field;
    size_t                *d_offsets;

    public:
      Field(size_t L, size_t T);
      Field(Element const &value, size_t L, size_t T);
      Field(Field< Element > const &other);

      template< typename Super >
      Field(Component< Super, Element> const &component); // Creation of a field through a Component of another field.

      template< typename Super >
      Field(hcComponent< Super, Element> const &component); // Creation of a field through a hcComponent of another field.

      explicit Field(hcField< Element  > const &other);
      Field< Element > &operator=(Field< Element > const &other);

      ~Field();

#include "Field/Field.iterator"
#include "Field/Field.const_iterator"

      size_t L() const;
      size_t T() const;
      size_t volume() const;

      double weave(Base::weaveOperator wea_OP, double nodeval) const;
      std::complex < double > weave(Base::weaveOperator wea_OP, std::complex < double > nodeval) const;

      iterator begin();
      iterator end();

      const_iterator begin() const;
      const_iterator end() const;

      template< typename Atom >
      Component< Element, Atom > component(size_t const component);

      Field< Element > &shift(Base::SpaceTimeIndex const idx, Base::Direction const shift);

#include "Field/Field.operators"

      hcField< Element > dagger() const;

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

      size_t size() const;
      size_t spatialSize() const; // (for looping over a single timeslice)
      void isolate();
      void refCountUp();
      void destroy();
      void fill(Element const &element); // Flush a field with a constant quantity

    private:
      size_t shiftIdxToZero(size_t const idx) const;
      size_t shiftIdxToOffset(size_t const idx) const;
  };

  template< typename Element >
  class hcField
  {
    friend class Field< Element >;
    Field< Element > const &d_parent;

    hcField(Field< Element > const &parent);

    public:
      Field< Element > const &dagger() const;

      typename Field< Element >::const_iterator begin() const;
      typename Field< Element >::const_iterator end() const;

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
