#ifndef GUARD_CORE_COMPONENT_H
#define GUARD_CORE_COMPONENT_H

namespace Core
{
  template< typename Element, size_t L, size_t T, typename Atom >
  class hcComponent;

  template< typename Element, size_t L, size_t T, typename Atom >
  class Component
  {
    friend class Field< Element, L, T >;

    Field< Element, L, T > &d_parent;
    SpaceTimeIndex const d_component;

    public:
      Component(Field< Element, L, T > &parent, Core::SpaceTimeIndex component);

#include "Component.iterator"

      iterator begin();
      iterator end();

      const_iterator begin() const;
      const_iterator end() const;

      Atom &element(Core::SpaceTimeIndex const *idx);
      hcComponent< Element, L, T, Atom > const dagger() const;

      #include "Component.operators"
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

#include "Component.inlines"
#include "Component.iterator.inlines"
#include "hcComponent.inlines"

#endif
