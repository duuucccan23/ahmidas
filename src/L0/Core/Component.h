#pragma once

#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>

namespace Core
{
  template< typename Element, typename Atom >
  class hcComponent;

  template< typename Element, typename Atom >
  class Component
  {
    friend class Field< Element >;

    Field< Element > &d_parent;
    size_t const      d_component;

    public:
      Component(Field< Element > &parent, size_t const component);

      #include "Component/Component.iterator"

      iterator begin();
      iterator end();

      size_t L() const;
      size_t T() const;
      size_t volume() const;

      const_iterator begin() const;
      const_iterator end() const;

      hcComponent< Element, Atom > const dagger() const;

      Atom &getMemoryIndex(size_t const idx);
      Atom const &getMemoryIndex(size_t const idx) const;

      #include "Component/Component.operators"

    private:
      Atom &getPhysicalIndex(size_t const idx);
      Atom const &getPhysicalIndex(size_t const idx) const;
  };

  template< typename Element, typename Atom >
  class hcComponent
  {
    friend class Component< Element, Atom >;

    Component< Element, Atom > const &d_parent;


    public:
      hcComponent(Component< Element, Atom > const &parent);
      Component< Element, Atom > const &dagger() const; //NOTE We wanted to move this

      size_t L() const;
      size_t T() const;
      size_t volume() const;

      typename Component< Element, Atom >::const_iterator &begin() const;
      typename Component< Element, Atom >::const_iterator &end() const;
      Atom const &parentIdx(size_t const idx) const;

  };
}

#include "Component/Component.inlines"
#include "Component/Component.iterator.inlines"
#include "Component/hcComponent.inlines"
