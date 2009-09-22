#pragma once

#include <L0/Base/Base.h>

namespace Source
{
  template< size_t L, size_t T, Base::SourceType Type >
  class Point
  {};

  template< size_t L, size_t T >
  class Point< L, T, (sou_SINGLET | sou_UNPOLARIZED) >
  {
    size_t d_coord[3];

    public:
      Point(size_t const x, size_t const y, size_t const z);
      Point(size_t const *coord);

      size_t coord(Base::SpaceTimeIndex const idx) const;
      size_t const *coord() const;
  };

  template< size_t L, size_t T >
  class Point< L, T, (sou_SINGLET | sou_PARTLY_POLARIZED) >
  {
    size_t                 d_coord[3];
    std::complex< double > d_field[4];

    public:
      Point(size_t const x, size_t const y, size_t const z);
      Point(size_t const *coord);

      size_t coord(Base::SpaceTimeIndex const idx) const;
      size_t const *coord() const;
      
      std::complex< double > *field();
      std::complex< double > const *field() const;
  };

  template< size_t L, size_t T >
  class Point< L, T, (sou_SINGLET | sou_FULLY_POLARIZED) >
  {
    size_t           d_coord[4];
    Base::DiracIndex d_index;

    public:
      Point(size_t const x, size_t const y, size_t const z);
      Point(size_t const *coord);
      size_t coord(Base::SpaceTimeIndex const idx) const;
      size_t const *coord() const;
      
      Base::DiracIndex &direction();
      Base::DiracIndex const &direction() const;
  };
  
  template< size_t L, size_t T >
  class Point< L, T, (sou_TRIPLET | sou_UNPOLARIZED) >
  {
    size_t                 d_coord[4];
    std::complex< double > d_field[3];

    public:
      Point(size_t const x, size_t const y, size_t const z);
      Point(size_t const *coord);

      size_t coord(Base::SpaceTimeIndex const idx) const;
      size_t const *coord() const;
      
      std::complex< double > *field();
      std::complex< double > const *field() const;
  };

  template< size_t L, size_t T >
  class Point< L, T, (sou_TRIPLET | sou_PARTLY_POLARIZED) >
  {
    size_t                 d_coord[4];
    std::complex< double > d_field[12];

    public:
      Point(size_t const x, size_t const y, size_t const z);
      Point(size_t const *coord);

      size_t coord(Base::SpaceTimeIndex const idx) const;
      size_t const *coord() const;
      
      std::complex< double > *field();
      std::complex< double > const *field() const;
  };

  template< size_t L, size_t T >
  class Point< L, T, (sou_TRIPLET | sou_FULLY_POLARIZED) >
  {
    size_t           d_coord[4];
    Base::DiracIndex d_index;
    std::complex< double > d_field[3];

    public:
      Point(size_t const x, size_t const y, size_t const z);
      Point(size_t const *coord);
      size_t coord(Base::SpaceTimeIndex const idx) const;
      size_t const *coord() const;
      
      std::complex< double > *field();
      std::complex< double > const *field() const;
  };
}

#include "Point/Point.inlines"
#include "Point/Point_Point.template"
