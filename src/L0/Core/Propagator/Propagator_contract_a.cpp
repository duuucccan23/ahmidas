#include "Propagator.ih"

namespace Core
{
  Field< Dirac::Matrix > *Propagator::contract(Propagator const &other) const
  {
    assert (T()==other.T() && L()==other.L());

    Field< Dirac::Matrix > *field = new Field< Dirac::Matrix > (L(), T());

    Propagator::const_iterator Ia(begin());
    Propagator::const_iterator Ib(other.begin());

    Field< Dirac::Matrix >::iterator Ic(field->begin());

    Dirac::Matrix RR, GG, BB;

    while(Ia != end())
    {
      QCD::Tensor tmp(*Ia);
      tmp.leftMultiply(*Ib);
      tmp.getDiracMatrix(RR, Base::col_RED,   Base::col_RED);
      tmp.getDiracMatrix(GG, Base::col_GREEN, Base::col_GREEN);
      tmp.getDiracMatrix(BB, Base::col_BLUE,  Base::col_BLUE);
      // *Ic should be initialized to zero by standard constructor of Dirac::Matrix
      (*Ic) += RR;
      (*Ic) += GG;
      (*Ic) += BB;
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
