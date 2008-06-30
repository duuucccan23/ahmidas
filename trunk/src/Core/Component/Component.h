#ifndef GUARD_CORE_COMPONENT_H
#define GUARD_CORE_COMPONENT_H

namespace Core
{
  template< typename Element, size_t L, size_t T, typename Atom >
  class hcComponent;
  
  template< typename Element, size_t L, size_t T, typename Atom >
  class Component
  {
    Field< Element, L, T > &d_parent;
    Core::SpaceTimeIndex    d_component;

    Component(Field< Element, L, T > &parent, Core::SpaceTimeIndex component);

#include "Component.iterator"

    public:
      iterator begin();
      iterator end();

      Atom &element(Core::SpaceTimeIndex const *idx);
      hcComponent< Element, L, T, Atom > const &dagger() const;
      
      #include "Component.operators"
  };
  
  template< typename Element, size_t L, size_t T, typename Atom >
  class hcComponent
  {
    Component< Element, L, T, Atom > const &d_parent;
    
    hcComponent(Component< Element,  L,  T,  Atom > const &parent);
    
    public:
      Component< Element, L, T, Atom > const &dagger() const;
  };
  
}

#include "Component.inlines"
#include "Component.iterator.inlines"
#include "hcComponent.inlines"

#endif
