#include "APE.ih"

namespace Smear
{
  void APE::smear(Core::Field< QCD::Gauge > &field, size_t const iterations) const
  {
    for (size_t ctr = 0; ctr < iterations; ctr++)
      smear(field);
  }

  void APE::smear(Core::Field< QCD::Gauge > &field, size_t const iterations, size_t const timeslice) const
  {
    size_t const L = field.L();
    size_t const T = field.T();

    assert(timeslice < T);
    Base::Weave weaveAll(L, T);
    assert(weaveAll.rank() == 0);
    weaveAll.barrier();
    Base::Weave weaveTimeslice(L, 1);

    Core::Field< QCD::Gauge > timeslice_field(L, 1);

    // copy one timeslice of the field and smear that

    for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L; idx_X++)
        {
          timeslice_field[weaveTimeslice.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, 1)] =
            field[weaveAll.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, timeslice)];
        }
      }
    }


    for (size_t ctr = 0; ctr < iterations; ctr++)
      smear(timeslice_field);


    for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L; idx_X++)
        {
          field[weaveAll.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, timeslice)]
           = timeslice_field[weaveTimeslice.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, 1)];
        }
      }
    }
  }

}
