#ifndef GUARD_CORE_H
#define GUARD_CORE_H

#include <algorithm>
#include <cassert>
#include <stdint.h>

namespace Core
{
  // Definitions for clarity in refering to indices -- beware the inverted order...
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
    col_RED = 0,
    col_GREEN = 1,
    col_BLUE = 2
  };

  enum DiracIndex
  {
    gam_1 = 0,
    gam_2 = 1,
    gam_3 = 2,
    gam_4 = 3
  };

  enum MpiTags
  {
    TAG_GAUGEFIELD,
    TAG_FILE_DISTRIBUTION
  };
}

#include "Core/Core.inlines"

#endif
