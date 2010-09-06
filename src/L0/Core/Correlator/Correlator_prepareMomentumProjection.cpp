#include "Correlator.ih"

namespace Core
{
  // performs a summation over timeslices including non-zero momentum projection
  void Correlator::prepareMomentumProjection(int const * const position_offset)
  {
    size_t const L(this->L());
    size_t const T(this->T());

    if (s_xRelative == NULL)
      s_xRelative =  new Field< int * >(L , T);
    else if (s_xRelative->L() != L ||  s_xRelative->T() != T)
    {
      for (Field< int *>::iterator I(s_xRelative->begin()); I != s_xRelative->end(); ++I)
        delete [] (*I);
      delete s_xRelative; // we can do this here because it is a static pointer,
      // this will affect all instances of Core::Correlator
      s_xRelative =  new Field< int * >(L , T);
    }

    size_t localIndex;
    for(size_t idx_T = 0; idx_T < T; idx_T++)
    {
      for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
      {
        for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
        {
          for(size_t idx_X = 0; idx_X < L; idx_X++)
          {
            localIndex = d_weave->globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);

            if (localIndex == d_weave->localVolume())
              continue;

            (*s_xRelative)[localIndex] = new int[3];
            ((*s_xRelative)[localIndex])[Base::idx_X] = (int(idx_X) - position_offset[Base::idx_X]) % T;
            ((*s_xRelative)[localIndex])[Base::idx_Y] = (int(idx_Y) - position_offset[Base::idx_Y]) % T;
            ((*s_xRelative)[localIndex])[Base::idx_Z] = (int(idx_Z) - position_offset[Base::idx_Z]) % T;
          }
        }
      }
    }
  }
}