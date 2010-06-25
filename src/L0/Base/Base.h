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
    sou_UNPOLARIZED      = 0,
    sou_PARTLY_POLARIZED = 1,
    sou_FULLY_POLARIZED  = 2
  };

  enum SourceColorState
  {
    sou_WHITE   = 0,
    sou_GENERIC = 1,
    sou_PURE    = 2
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
    proj_PARITY_PLUS_TM_STAR,
    proj_PARITY_MINUS_TM,
    proj_NO_PROJECTOR // this one does not do anything
  };

  enum Operator
  {
    op_UNITY = -1,
    op_GAMMA_4 = 0,
    op_GAMMA_1 = 1,
    op_GAMMA_2 = 2,
    op_GAMMA_3 = 3,
    op_GAMMA_5 = 4,
    op_GAMMA_4_GAMMA_5 = 5,
    op_GAMMA_1_GAMMA_5 = 6,
    op_GAMMA_2_GAMMA_5 = 7,
    op_GAMMA_3_GAMMA_5 = 8,
    op_CONSERVED_GAMMA_4 = 16,
    op_O44 = 32,
    op_O11 = 33,
    op_O22 = 34,
    op_O33 = 35
  };

  extern bool const bigEndian;
}

#include "Base/Base.inlines"
