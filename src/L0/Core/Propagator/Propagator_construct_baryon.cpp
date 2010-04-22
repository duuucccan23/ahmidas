#include "Propagator.ih"

namespace Core
{

  Field< Dirac::Matrix > *Propagator::construct_baryon(Propagator const &no2, Propagator const &no3,
                                                       Base::BaryonInterpolatingField const ipol) const
  {
    assert (T()==no2.T() && L()==no2.L() && T()==no3.T() && L()==no3.L());

    Field< Dirac::Matrix > *field = new Field< Dirac::Matrix > (L(), T());

    Propagator::const_iterator Ia(begin());
    Propagator::const_iterator Ib(no2.begin());
    Propagator::const_iterator Ic(no3.begin());

    Field< Dirac::Matrix >::iterator Id(field->begin());

    while(Ia != this->end())
    {
//        std::cout << (*Ia) << std::endl;

      // the constructor of reducedTensor has to process
      // the (anti-symmetric) colour structure epsilon_abc
      // as well as the interpolating Dirac structure
      // between the second and third quark field
      QCD::getDiracMatrix((*Id),(*Ia),(*Ib),(*Ic), ipol);
      ++Ia;
      ++Ib;
      ++Ic;
      ++Id;
    }
    assert(Ia==end());
    assert(Ib==no2.end());
    assert(Ic==no3.end());
    assert(Id==field->end());
    return field;
  }

}
