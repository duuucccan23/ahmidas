#ifndef GUARD_CORE_COMPONENT_H
#define GUARD_CORE_COMPONENT_H

namespace Core
{
  template< typename Element, size_t L, size_t T, typename Atom >
  class Component
  {
    Field< Element, L, T > &d_parent;
    short                   d_component;

    Component(Field< Element, L, T > &parent, short component);

#include "Component.iterator"

    public:
      iterator begin();
      iterator end();

      Atom &element(short const *idx);

      template< typename Scalar >
      void leftMultiply(Scalar const &other);

      template< typename OtherElement >
      void leftMultiply(Field< OtherElement, L, T > const &other);

      template< typename OtherElement, typename OtherAtom  >
      void leftMultiply(Component< OtherElement, L, T, OtherAtom > const &other);

      template< typename Scalar >
      void rightMultiply(Scalar const &other);

      template< typename OtherElement >
      void rightMultiply(Field< OtherElement, L, T > const &other);

      template< typename OtherElement, typename OtherAtom  >
      void rightMultiply(Component< OtherElement, L, T, OtherAtom > const &other);

      template< typename OtherElement >
      Component< Element, L, T, Atom > &operator+=(Field< OtherElement, L, T > const &other);

      template< typename OtherElement, typename OtherAtom  >
      Component< Element, L, T, Atom > &operator+=(Component< OtherElement, L, T, OtherAtom > const &other);

      template< typename OtherElement >
      Component< Element, L, T, Atom > &operator-=(Field< OtherElement, L, T > const &other);

      template< typename OtherElement, typename OtherAtom  >
      Component< Element, L, T, Atom > &operator-=(Component< OtherElement, L, T, OtherAtom > const &other);

      template< typename Scalar >
      Component< Element, L, T, Atom > &operator*=(Scalar const &rhand);

      template< typename OtherElement >
      Component< Element, L, T, Atom > &operator*=(Field< OtherElement, L, T > const &field);

      template< typename OtherElement, typename OtherAtom >
      Component< Element, L, T, Atom > &operator*=(Component< OtherElement, L, T, OtherAtom > const &field);

      template< typename Scalar >
      Component< Element, L, T, Atom > &operator/=(Scalar const &rhand);

      template< typename OtherElement >
      Component< Element, L, T, Atom > &operator/=(Field< OtherElement, L, T > const &field);

      template< typename OtherElement, typename OtherAtom >
      Component< Element, L, T, Atom > &operator/=(Component< OtherElement, L, T, OtherAtom > const &field);
  };
}

#include "Component.inlines"
#include "Component.iterator.inlines"

#endif
