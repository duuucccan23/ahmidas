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
    sou_UNPOLARIZED = 0,
    sou_PARTLY_POLARIZED = 1,
    sou_FULLY_POLARIZED = 2
  };

  enum SourceColorState
  {
    sou_WHITE = 0,
    sou_PURE = 1,
    sou_GENERIC = 2
  };

  extern bool const bigEndian;
}

#include "Base/Base.inlines"
