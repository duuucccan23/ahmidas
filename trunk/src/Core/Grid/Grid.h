#ifndef GUARD_CORE_GRID_H
#define GUARD_CORE_GRID_H

#include <algorithm>
#include <cerrno>
#include <functional>
#include <iostream>

#include <mpi/mpi.h>

#include "../Core.h"

namespace Core
{
  template< size_t L, size_t T >
  class Grid
  {
    MPI::Cartcomm d_grid;
    MPI::Cartcomm d_timeSlice;
    MPI::Cartcomm d_backbone;

    size_t d_dims[4]; //number of nodes in each dimension
    size_t d_sizes[4]; //local dimension sizes
    size_t d_coords[4];
    size_t d_surfaces[4];
    size_t d_dimSizes[4];

    size_t d_localVolume;
    size_t d_bufferVolume;
    size_t d_contiguousBlock;

    bool d_bigEndian;

    public:
      Grid();
      ~Grid();

      MPI::Cartcomm &grid();
      MPI::Cartcomm &timeSlice();
      MPI::Cartcomm &backbone();

      size_t rank() const; //rank of this node
      size_t rank(size_t index) const; //rank of the node holding the lattice site with index index
//      size_t rank(size_t const *site) const; //rank of the node holding the lattice site with coordinates site
      size_t rank(int const *coords) const; //rank of the node with node coordinates coords

      size_t const *coords() const; //node coordinates of this node
      size_t coord(size_t idx) const; //node coordinate in dimension idx

      size_t const *dims() const; //number of nodes in each dimension
      size_t dim(size_t idx) const; //number of nodes in dimension idx

      size_t const *sizes() const; //dimension sizes on this (=every) node
      size_t size(size_t idx) const; //dimension idx size on this (=every) node

      size_t const *surfaces() const;
      size_t surface(size_t idx) const;

      size_t localVolume() const;
      size_t totalVolume() const;

      size_t const *dimSizes() const;
      size_t dimSize(size_t idx) const;

      size_t neighbour(SpaceTimeIndex idx, Direction dir) const;

      size_t bufferVolume() const;
      size_t contiguousBlock() const;

      bool bigEndian() const;

    private:
      Grid(const Grid& other); // Sterile constructor
      void initContiguousBlock();
      size_t greatestCommonDivisor(size_t x, size_t y);
  };
}

#include "Grid.inlines"
#include "Grid_a.template"
#include "Grid_b.template"
#include "Grid_c.template"
#include "Grid_d.template"

#endif
