#include "Propagator.ih"

namespace Core
{

  void Propagator::gaugeTransform_fixedSource(Field< SU3::Matrix > const &gaugeTrafo, size_t const * sourcePos)
  {
    assert(gaugeTrafo.L() == L() && gaugeTrafo.T() == T());
    Base::Weave weave(L(), T());
    SU3::Matrix matrixAtSource;

    size_t rank = weave.rank(sourcePos);

    size_t localIndex = weave.globalCoordToLocalIndex(sourcePos[Base::idx_X],
                                                      sourcePos[Base::idx_Y],
                                                      sourcePos[Base::idx_Z],
                                                      sourcePos[Base::idx_T]);
    if (rank == weave.rank())
    {
      matrixAtSource = SU3::Matrix((gaugeTrafo[localIndex]).dagger());
    }

    weave.broadcast(&matrixAtSource, 1, rank);

    Propagator::iterator I_p = begin();
    Field< SU3::Matrix >::const_iterator I_g = gaugeTrafo.begin();
    while (I_p != end())
    {
      (*I_p).leftMultiply(matrixAtSource).rightMultiply(*I_g);
      ++I_p;
      ++I_g;
    }
    assert(I_g == gaugeTrafo.end());

  }

}