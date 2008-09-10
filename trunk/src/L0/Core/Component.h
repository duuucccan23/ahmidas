#ifndef GUARD_CORE_COMPONENT_H
#define GUARD_CORE_COMPONENT_H

#include <L0/Core/Field.h>

namespace Core
{
  template< typename Element, size_t L, size_t T, typename Atom >
  class hcComponent;

  template< typename Element, size_t L, size_t T, typename Atom >
  class Component
  {
    friend class Field< Element, L, T >;

    Unifold< L, T > &d_unifold;
    Field< Element, L, T > &d_parent;
    size_t const d_component;

    public:
      Component(Field< Element, L, T > &parent, size_t const component);

      #include "Component/Component.iterator"

      iterator begin();
      iterator end();

      const_iterator begin() const;
      const_iterator end() const;

      hcComponent< Element, L, T, Atom > const dagger() const;

      Atom &getMemoryIndex(size_t const idx);
      Atom const &getMemoryIndex(size_t const idx) const;

      #include "Component/Component.operators"

    private:
      Atom &getPhysicalIndex(size_t const idx);
      Atom const &getPhysicalIndex(size_t const idx) const;
  };

  template< typename Element, size_t L, size_t T, typename Atom >
  class hcComponent
  {
    friend class Component< Element, L, T, Atom >;

    Component< Element, L, T, Atom > const &d_parent;

    hcComponent(Component< Element,  L,  T,  Atom > const &parent);

    public:
      Component< Element, L, T, Atom > const &dagger() const;

      typename Component< Element, L, T, Atom >::const_iterator &begin() const;
      typename Component< Element, L, T, Atom >::const_iterator &end() const;
  };
}

#include "Component/Component.inlines"
#include "Component/Component.iterator.inlines"
#include "Component/hcComponent.inlines"

#endif
