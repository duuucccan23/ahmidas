#ifndef GUARD_CORE_COMPONENT_H
#define GUARD_CORE_COMPONENT_H

namespace Core
{
  template< typename Element, typename Atom >
  class Component
  {
    Field< Element > &d_parent;
    short             d_component;

    Component(Field< Element > &parent, short component);
    
#include "Component.iterator"
    
    public:
      iterator begin();
      iterator end();
      
      Atom &element(short const *idx);

      template< typename Scalar >
      void leftMultiply(Scalar const &other);
      
      template< typename OtherElement >
      void leftMultiply(Field< OtherElement > const &other);          
      
      template< typename OtherElement, typename OtherAtom  >
      void leftMultiply(Component< OtherElement, OtherAtom > const &other);
      
      template< typename Scalar >
      void rightMultiply(Scalar const &other);
      
      template< typename OtherElement >
      void rightMultiply(Field< OtherElement > const &other);          
      
      template< typename OtherElement, typename OtherAtom  >
      void rightMultiply(Component< OtherElement, OtherAtom > const &other);         
      
      template< typename OtherElement >
      Component< Element, Atom > &operator+=(Field< OtherElement > const &other);      
      
      template< typename OtherElement, typename OtherAtom  >
      Component< Element, Atom > &operator+=(Component< OtherElement, OtherAtom > const &other);
      
      template< typename OtherElement >
      Component< Element, Atom > &operator-=(Field< OtherElement > const &other);          
      
      template< typename OtherElement, typename OtherAtom  >
      Component< Element, Atom > &operator-=(Component< OtherElement, OtherAtom > const &other);
    
      template< typename Scalar >
      Component< Element, Atom > &operator*=(Scalar const &rhand);      

      template< typename OtherElement >
      Component< Element, Atom > &operator*=(Field< OtherElement > const &field); 
      
      template< typename OtherElement, typename OtherAtom >
      Component< Element, Atom > &operator*=(Component< OtherElement, OtherAtom > const &field);       
      
      template< typename Scalar >
      Component< Element, Atom > &operator/=(Scalar const &rhand);

      template< typename OtherElement >
      Component< Element, Atom > &operator/=(Field< OtherElement > const &field); 
      
      template< typename OtherElement, typename OtherAtom >
      Component< Element, Atom > &operator/=(Component< OtherElement, OtherAtom > const &field);      
  };
}

#include "Component.inlines"
#include "Component.iterator.inlines"
#include "Component.templates"

#endif
