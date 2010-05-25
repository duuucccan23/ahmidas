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


/*
    op_GAMMA_4 = 0
    op_GAMMA_1 = 1
    op_GAMMA_2 = 2
    op_GAMMA_3 = 3

    op_O44 = 32
    op_O11 = 33
    op_O22 = 34
    op_O33 = 35*/


  void Propagator::operator*=(Base::Operator const& O)
  {
    switch (O)
    {
      case Base::op_GAMMA_4:
      {
        Dirac::Gamma< 4 > gamma4;
        (*this) *= gamma4;
        break;
      }
      case Base::op_GAMMA_1:
      {
        Dirac::Gamma< 1 > gamma1;
        (*this) *= gamma1;
        break;
      }
      case Base::op_GAMMA_2:
      {
        Dirac::Gamma< 2 > gamma2;
        (*this) *= gamma2;
        break;
      }
      case Base::op_GAMMA_3:
      {
        Dirac::Gamma< 3 > gamma3;
        (*this) *= gamma3;
        break;
      }
//       case Base::op_O44:
//       {
//         Dirac::Gamma< 4 > gamma4;
//         (*this) *= gamma4;
//         break;
//       }
      default:
        std::cerr << "Error in void Propagator::operator*=(Base::Operator const& O)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        exit(1);
    }
  }


// Propagator operator*(Base::Operator const& O, Propagator const& P)
// {
//   Propagator tmp(P);
// 
//   return tmp;
// }


}
