#pragma once

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdint.h>

namespace Base
{
  // Definitions for clarity in refering to indices
  enum SpaceTimeIndex
  {
    idx_X = 0,
    idx_Y = 1,
    idx_Z = 2,
    idx_T = 3
  };

  // Definitions to be used in refering to directions for shifts
  enum Direction
  {
    dir_DOWN = -1,
    dir_UP = 1
  };

  enum ColourIndex
  {
    col_RED   = 0,
    col_GREEN = 1,
    col_BLUE  = 2
  };

  enum DiracIndex
  {
    gam_1 = 0,
    gam_2 = 1,
    gam_3 = 2,
    gam_4 = 3
  };

  // Definitions to be used in the characterization of sources
  enum SourcePolarization
  {
    sou_UNPOLARIZED,
    sou_PARTLY_POLARIZED,
    sou_FULLY_POLARIZED
  };

  enum SourceColorState
  {
    sou_WHITE,
    sou_PURE,
    sou_GENERIC
  };

  enum SourceSpatialLayout
  {
    sou_POINT,
    sou_FULL,
    sou_WALL
  };

  enum weaveOperator
  {
    wea_SUM,
    wea_XOR
  };

//   enum boundaryConditions
//   {
//     bc_ANTIPERIODIC_FIXED;
//     bc_ANTIPERIODIC_UNIFORM;
//     bc_PERIODIC_FIXED;
//     bc_PERIODIC_UNIFORM;
//   }

  // Definitions to be used for baryon contractions
  enum BaryonInterpolatingField
  {
    bar_PROTON
  };

  enum BaryonPropagatorProjector
  {
    proj_PARITY_PLUS_STD,
    proj_PARITY_MINUS_STD,
    proj_PARITY_PLUS_TM,
    proj_PARITY_MINUS_TM
  };

  extern bool const bigEndian;
}

#include "Base/Base.inlines"
