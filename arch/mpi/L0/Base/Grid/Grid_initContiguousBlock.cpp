#include "Grid.ih"

namespace Base
{
  void Grid::initContiguousBlock()
  {
    size_t cutIdx = 0;
    while ((d_dims[cutIdx] == 1) && (cutIdx < 4))
      ++cutIdx;
    d_contiguousBlock = dimSize(cutIdx);
  }
}