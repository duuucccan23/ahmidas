#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

namespace Core
{

  // for 1-derivative operators or conserved vector current
  Field< Dirac::Matrix > Propagator::contractWithOperatorInsertion(Base::Operator const O, Field< QCD::Gauge > * const gauge_field, Propagator const &fromRight, Core::Field< Dirac::Matrix > &timesG5)
  {

    assert(gauge_field != NULL);

    Core::Field< Dirac::Matrix > result(L(),T());

    Dirac::Gamma< 5 > gamma5;

    switch (O)
    {
      case Base::op_CONSERVED_GAMMA_4:
      {
        Dirac::Gamma< 4 > gamma4;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        Propagator copy(*this);
        // gamma_mu - identity
        copy *= gamma4;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  -= (*this);
        copy5 -= (*this);
        result = copy.contract(tmp);
        timesG5 = copy5.contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator copy(*this);
        // gamma_mu + identity
        copy *= gamma4;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  += (*this);
        copy5 += (*this);
        (*result) += (*(copy.contract(tmp)));
        (*timesG5) += (*(copy5.contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma4;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp += (*this);
        tmp5 += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp5.shift(Base::idx_T, Base::dir_DOWN);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        (*result) += (*(tmp.contract(fromRight)));
        (*timesG5) += (*(tmp5.contract(fromRight)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma4;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp -= (*this);
        tmp5 -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp5.shift(Base::idx_T, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        (*timesG5) += *(tmp5.contract(fromRight));
        }
        break;
      }
      case Base::op_CONSERVED_GAMMA_1:
      {
        Dirac::Gamma< 1 > gamma1;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        Propagator copy(*this);
        // gamma_mu - identity
        copy *= gamma1;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  -= (*this);
        copy5 -= (*this);
        result = copy.contract(tmp);
        timesG5 = copy5.contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator copy(*this);
        // gamma_mu + identity
        copy *= gamma1;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  += (*this);
        copy5 += (*this);
        (*result) += (*(copy.contract(tmp)));
        (*timesG5) += (*(copy5.contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma1;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp += (*this);
        tmp5 += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp5.shift(Base::idx_X, Base::dir_DOWN);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        (*result) += (*(tmp.contract(fromRight)));
        (*timesG5) += (*(tmp5.contract(fromRight)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma1;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp -= (*this);
        tmp5 -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp5.shift(Base::idx_X, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        (*timesG5) += *(tmp5.contract(fromRight));
        }
        break;
      }
      case Base::op_CONSERVED_GAMMA_2:
      {
        Dirac::Gamma< 2 > gamma2;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        Propagator copy(*this);
        // gamma_mu - identity
        copy *= gamma2;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  -= (*this);
        copy5 -= (*this);
        result = copy.contract(tmp);
        timesG5 = copy5.contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator copy(*this);
        // gamma_mu + identity
        copy *= gamma2;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  += (*this);
        copy5 += (*this);
        (*result) += (*(copy.contract(tmp)));
        (*timesG5) += (*(copy5.contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma2;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp += (*this);
        tmp5 += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp5.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        (*result) += (*(tmp.contract(fromRight)));
        (*timesG5) += (*(tmp5.contract(fromRight)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma2;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp -= (*this);
        tmp5 -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp5.shift(Base::idx_Y, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        (*timesG5) += *(tmp5.contract(fromRight));
        }
        break;
      }
      case Base::op_CONSERVED_GAMMA_3:
      {
        Dirac::Gamma< 3 > gamma3;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        Propagator copy(*this);
        // gamma_mu - identity
        copy *= gamma3;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  -= (*this);
        copy5 -= (*this);
        result = copy.contract(tmp);
        timesG5 = copy5.contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator copy(*this);
        // gamma_mu + identity
        copy *= gamma3;
        Propagator copy5(copy);
        copy5 *= gamma5;
        copy  += (*this);
        copy5 += (*this);
        (*result) += (*(copy.contract(tmp)));
        (*timesG5) += (*(copy5.contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma3;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp += (*this);
        tmp5 += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp5.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        (*result) += (*(tmp.contract(fromRight)));
        (*timesG5) += (*(tmp5.contract(fromRight)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma3;
        Propagator tmp5(tmp);
        tmp5 *= gamma5;
        tmp -= (*this);
        tmp5 -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        (tmp5.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp5.shift(Base::idx_Z, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        (*timesG5) += *(tmp5.contract(fromRight));
        }
        break;
      }

      // 1-derivative operators

      case Base::op_O44:
      {
        Dirac::Gamma< 4 > gamma4;
        (*this) *= gamma4;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O11:
      {
        Dirac::Gamma< 1 > gamma1;
        (*this) *= gamma1;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O22:
      {
        Dirac::Gamma< 2 > gamma2;
        (*this) *= gamma2;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O33:
      {
        Dirac::Gamma< 3 > gamma3;
        (*this) *= gamma3;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O41:
      {
        Dirac::Gamma< 4 > gamma4;
        (*this) *= gamma4;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O42:
      {
        Dirac::Gamma< 4 > gamma4;
        (*this) *= gamma4;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O43:
      {
        Dirac::Gamma< 4 > gamma4;
        (*this) *= gamma4;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O14:
      {
        Dirac::Gamma< 1 > gamma1;
        (*this) *= gamma1;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O12:
      {
        Dirac::Gamma< 1 > gamma1;
        (*this) *= gamma1;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O13:
      {
        Dirac::Gamma< 1 > gamma1;
        (*this) *= gamma1;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O24:
      {
        Dirac::Gamma< 2 > gamma2;
        (*this) *= gamma2;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O21:
      {
        Dirac::Gamma< 2 > gamma2;
        (*this) *= gamma2;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O23:
      {
        Dirac::Gamma< 2 > gamma2;
        (*this) *= gamma2;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O34:
      {
        Dirac::Gamma< 3 > gamma3;
        (*this) *= gamma3;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O31:
      {
        Dirac::Gamma< 3 > gamma3;
        (*this) *= gamma3;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O32:
      {
        Dirac::Gamma< 3 > gamma3;
        (*this) *= gamma3;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }

      default:
        std::cerr << "Error in void Propagator::contractWithOperatorInsertion(...)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        delete result;
        delete timesG5;
        exit(1);
    }
    return result;


  }
}


/* prototype for O_mu_nu

      case Base::op_OMUNU:
      {
        Dirac::Gamma< MU > gammaMU;
        (*this) *= gammaMU;
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_NU, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_NU));
        Propagator tmp5(gamma5*tmp);
        timesG5 = contract(tmp5);
        result = contract(tmp);
        }
        {
        Propagator tmp(fromRight);
        tmp.isolate();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_NU).dagger());
        tmp.shift(Base::idx_NU, Base::dir_UP);
        Propagator tmp5(gamma5*tmp);
        (*timesG5) = (*contract(tmp5));
        (*result) -= (*(contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_NU, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_NU).dagger());
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) -= *(tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_NU));
        tmp.shift(Base::idx_NU, Base::dir_UP);
        Propagator tmp5(tmp*gamma5);
        (*timesG5) -= *(tmp5.contract(fromRight));
        (*result) += *(tmp.contract(fromRight));
        }
        break;
      }
*/








