#include "Propagator.ih"

namespace Core
{
  QCD::Tensor Propagator::operator()(size_t const * sinkSite) const
  {
    Base::Weave weave(L(), T());
    QCD::Tensor tensor;

    size_t rank = weave.rank(sinkSite);

    if (rank == weave.rank())
    {
      size_t localIndex = weave.globalCoordToLocalIndex(sinkSite[Base::idx_X],
                                                        sinkSite[Base::idx_Y],
                                                        sinkSite[Base::idx_Z],
                                                        sinkSite[Base::idx_T]);
      tensor = QCD::Tensor((*d_components)[localIndex]);
    }

    weave.broadcast(&tensor, 1, rank);

    return tensor;
  }


  void Propagator::operator*=(std::complex< double > const &factor)
  {
    isolate();
    Propagator::iterator I(begin());

    while(I != end())
    {
      (*I) *= factor;
      ++I;
    }
  }

  Field< Dirac::Matrix > *Propagator::operator*(Propagator const &other) const
  {
    assert (T()==other.T() && L()==other.L());

    Field< Dirac::Matrix > *field = new Field< Dirac::Matrix > (L(), T());

    Propagator::const_iterator Ia(begin());
    Propagator::const_iterator Ib(other.begin());

    Field< Dirac::Matrix >::iterator Ic(field->begin());

    while(Ia != end())
    {
      QCD::getDiracMatrix((*Ic), (*Ia),(*Ib));
      ++Ia;
      ++Ib;
      ++Ic;
    }
    assert(Ia==end());
    assert(Ib==other.end());
    assert(Ic==field->end());
    return field;
  }

}
