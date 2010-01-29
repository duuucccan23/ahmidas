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

  enum weaveOperator
  {
    wea_SUM,
    wea_XOR
  };

  extern bool const bigEndian;
}

#include "Base/Base.inlines"
