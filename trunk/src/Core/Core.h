#ifndef GUARD_CORE_H
#define GUARD_CORE_H

#include <algorithm>
#include <cassert>

namespace Core
{
  // Definitions for clarity in refering to indices -- beware the inverted order...
  enum SpaceTimeIndex
  {
    idx_Z = 0,
    idx_Y = 1,
    idx_X = 2,
    idx_T = 3
  };

  // Definitions to be used in refering to directions for shifts
  enum Direction
  {
    dir_DOWN = -1,
    dir_UP = 1
  };

  enum MpiTags
  {
    TAG_GAUGEFIELD
  };
}

#include "Core.inlines"

#endif
