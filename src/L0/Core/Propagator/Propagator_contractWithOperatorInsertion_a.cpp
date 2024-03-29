#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

namespace Core
{

  Core::Field< Dirac::Matrix > Propagator::contractWithOperatorInsertion(Base::Operator const O, Field< QCD::Gauge > const * const gauge_field, Propagator const &fromRight)
  {
    Core::Field< Dirac::Matrix > result(L(),T());
    switch (O)
    {
      case Base::op_GAMMA_4:
      {
        Dirac::Gamma< 4 > gamma4;
        (*this) *= gamma4;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_1:
      {
        Dirac::Gamma< 1 > gamma1;
        (*this) *= gamma1;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_2:
      {

        Dirac::Gamma< 2 > gamma2;
        (*this) *= gamma2;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_3:
      {
        Dirac::Gamma< 3 > gamma3;
        (*this) *= gamma3;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_45:
      {
        Dirac::Gamma< 45 > gamma4gamma5;
        (*this) *= gamma4gamma5;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_15:
      {
        Dirac::Gamma< 15 > gamma1gamma5;
        (*this) *= gamma1gamma5;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_25:
      {

        Dirac::Gamma< 25 > gamma2gamma5;
        (*this) *= gamma2gamma5;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_35:
      {
        Dirac::Gamma< 35 > gamma3gamma5;
        (*this) *= gamma3gamma5;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_12:
      {
        Dirac::Gamma< 12 > gamma1gamma2;
        (*this) *= gamma1gamma2;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_13:
      {
        Dirac::Gamma< 13 > gamma1gamma3;
        (*this) *= gamma1gamma3;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_14:
      {
        Dirac::Gamma< 14 > gamma1gamma4;
        (*this) *= gamma1gamma4;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_23:
      {
        Dirac::Gamma< 23 > gamma2gamma3;
        (*this) *= gamma2gamma3;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_24:
      {
        Dirac::Gamma< 24 > gamma2gamma4;
        (*this) *= gamma2gamma4;
        result = contract(fromRight);
        break;
      }
      case Base::op_GAMMA_34:
      {
        Dirac::Gamma< 34 > gamma3gamma4;
        (*this) *= gamma3gamma4;
        result = contract(fromRight);
        break;
      }
      case Base::op_CONSERVED_GAMMA_4:
      {
        assert(gauge_field != NULL);
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
        copy -= (*this);
        result = copy.contract(tmp);
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
        copy += (*this);
        result += (copy.contract(tmp));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma4;
        tmp += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        result += (tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma4;
        tmp -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        result += (tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_CONSERVED_GAMMA_1:
      {
        assert(gauge_field != NULL);
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
        copy -= (*this);
        result = copy.contract(tmp);
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
        copy += (*this);
        result += ((copy.contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma1;
        tmp += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        result += (tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma1;
        tmp -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        result += (tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_CONSERVED_GAMMA_2:
      {
        assert(gauge_field != NULL);
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
        copy -= (*this);
        result = copy.contract(tmp);
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
        copy += (*this);
        result += ((copy.contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma2;
        tmp += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        result += (tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma2;
        tmp -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        result += (tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_CONSERVED_GAMMA_3:
      {
        assert(gauge_field != NULL);
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
        copy -= (*this);
        result = copy.contract(tmp);
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
        copy += (*this);
        result += ((copy.contract(tmp)));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu + identity
        tmp *= gamma3;
        tmp += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        result += (tmp.contract(fromRight));
        }
        {
        Propagator tmp(*this);
        tmp.isolate();
        // gamma_mu - identity
        tmp *= gamma3;
        tmp -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        result += (tmp.contract(fromRight));
        }
        break;
      }
      case Base::op_O44:
      {

        //std::cout << "calculating product with operator O44 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 4 > gamma4;

        (*this) *= gamma4;

        {
        Propagator tmp(fromRight);
        //std::cout << "part 1 ... ";
        //std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        result = contract(tmp);
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        //std::cout << "part 2 ... ";
        //std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp.shift(Base::idx_T, Base::dir_UP);
        result -= ((contract(tmp)));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 3 ... ";
        //std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        result -= (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 4 ... ";
        //std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        result += (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }
        break;
      }
      case Base::op_O11:
      {
        //std::cout << "calculating product with operator O11 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 1 > gamma1;

        (*this) *= gamma1;

        {
        Propagator tmp(fromRight);
        //std::cout << "part 1 ... ";
        //std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        result = contract(tmp);
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        //std::cout << "part 2 ... ";
        //std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp.shift(Base::idx_X, Base::dir_UP);
        result -= ((contract(tmp)));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 3 ... ";
        //std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        result -= (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 4 ... ";
        //std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        result += (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }
        break;
      }
      case Base::op_O22:
      {
        //std::cout << "calculating product with operator O22 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 2 > gamma2;

        (*this) *= gamma2;

        {
        Propagator tmp(fromRight);
        //std::cout << "part 1 ... ";
        //std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        result = contract(tmp);
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        //std::cout << "part 2 ... ";
        //std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp.shift(Base::idx_Y, Base::dir_UP);
        result -= ((contract(tmp)));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 3 ... ";
        //std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        result -= (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 4 ... ";
        //std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        result += (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }

        break;
      }
      case Base::op_O33:
      {
        //std::cout << "calculating product with operator O33 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 3 > gamma3;

        (*this) *= gamma3;

        {
        Propagator tmp(fromRight);
        //std::cout << "part 1 ... ";
        //std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        result = contract(tmp);
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        //std::cout << "part 2 ... ";
        //std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp.shift(Base::idx_Z, Base::dir_UP);
        result -= ((contract(tmp)));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 3 ... ";
        //std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        result -= (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        //std::cout << "part 4 ... ";
        //std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        result += (tmp.contract(fromRight));
        //std::cout << "done." << std::endl;
        }

        break;
      }
      default:
        std::cerr << "Error in void Propagator::contractWithOperatorInsertion(...)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        exit(1);
    }
    return result;
  }

}
