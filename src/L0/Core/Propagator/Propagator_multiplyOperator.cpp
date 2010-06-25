#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

namespace Core
{

  Core::Field< Dirac::Matrix > *Propagator::contractWithOperatorInsertion(Base::Operator const O, Field< QCD::Gauge > * const gauge_field, Propagator const &fromRight)
  {
    Core::Field< Dirac::Matrix > *result = NULL;
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
      case Base::op_CONSERVED_GAMMA_4:
      {

        assert(gauge_field != NULL);
        // not implemented
        result = new Core::Field< Dirac::Matrix >(L(), T());
        // Dirac::Gamma< 4 > gamma4;
        // std::complex< double > ZERO(0, 0);
/*
        Propagator copy(*this);
        (*this) *= ZERO;

        {
        Propagator tmp(copy);
        std::cout.flush();
        tmp.isolate();
        // part 1 : times U_mu and shift up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        (*this) -= tmp;
        tmp *= gamma4;
        (*this) += tmp;
        }

        {
        Propagator tmp(copy);
        tmp.isolate();
        std::cout.flush();
        // part 2 : shift down and times U_mu^dagger
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        (*this) += tmp;
        tmp *= gamma4;
        (*this) += tmp;
        }
*/
        break;
      }
      case Base::op_O44:
      {

        std::cout << "calculating product with operator O44 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 4 > gamma4;

        (*this) *= gamma4;

        {
        Propagator tmp(fromRight);
        std::cout << "part 1 ... ";
        std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        result = contract(tmp);
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        std::cout << "part 2 ... ";
        std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        tmp.shift(Base::idx_T, Base::dir_UP);
        (*result) -= (*(contract(tmp)));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 3 ... ";
        std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        (*result) -= *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 4 ... ";
        std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }
        break;
      }
      case Base::op_O11:
      {
        std::cout << "calculating product with operator O11 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 1 > gamma1;

        (*this) *= gamma1;

        {
        Propagator tmp(fromRight);
        std::cout << "part 1 ... ";
        std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        result = contract(tmp);
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        std::cout << "part 2 ... ";
        std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        tmp.shift(Base::idx_X, Base::dir_UP);
        (*result) -= (*(contract(tmp)));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 3 ... ";
        std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        (*result) -= *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 4 ... ";
        std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }
        break;
      }
      case Base::op_O22:
      {
        std::cout << "calculating product with operator O22 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 2 > gamma2;

        (*this) *= gamma2;

        {
        Propagator tmp(fromRight);
        std::cout << "part 1 ... ";
        std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        result = contract(tmp);
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        std::cout << "part 2 ... ";
        std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        tmp.shift(Base::idx_Y, Base::dir_UP);
        (*result) -= (*(contract(tmp)));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 3 ... ";
        std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        (*result) -= *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 4 ... ";
        std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }

        break;
      }
      case Base::op_O33:
      {
        std::cout << "calculating product with operator O33 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 3 > gamma3;

        (*this) *= gamma3;

        {
        Propagator tmp(fromRight);
        std::cout << "part 1 ... ";
        std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        result = contract(tmp);
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(fromRight);
        tmp.isolate();
        std::cout << "part 2 ... ";
        std::cout.flush();
        // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        tmp.shift(Base::idx_Z, Base::dir_UP);
        (*result) -= (*(contract(tmp)));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 3 ... ";
        std::cout.flush();
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        (*result) -= *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 4 ... ";
        std::cout.flush();
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }

        break;
      }
      default:
        std::cerr << "Error in void Propagator::contractWithOperatorInsertion(...)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        delete result;
        exit(1);
    }
    return result;
  }


  void Propagator::leftMultiplyOperator(Base::Operator const O)
  {
    isolate();
    switch (O)
    {
      case Base::op_UNITY:
        break;
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
      case Base::op_GAMMA_4_GAMMA_5:
      {
        Dirac::Gamma< 45 > gamma4gamma5;
        (*this) *= gamma4gamma5;
        break;
      }
      default:
        std::cerr << "Error in void Propagator::multiplyOperator(Base::Operator const& O)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        exit(1);
    }
  }

  void Propagator::rightMultiplyOperator(Base::Operator const O)
  {
    isolate();
    switch (O)
    {
      case Base::op_UNITY:
        break;
      case Base::op_GAMMA_4:
      {
        Dirac::Gamma< 4 > gamma4;
        (*this).rightMultiply(gamma4);
        break;
      }
      case Base::op_GAMMA_1:
      {
        Dirac::Gamma< 1 > gamma1;
        (*this).rightMultiply(gamma1);
        break;
      }
      case Base::op_GAMMA_2:
      {

        Dirac::Gamma< 2 > gamma2;
        (*this).rightMultiply(gamma2);
        break;
      }
      case Base::op_GAMMA_3:
      {
        Dirac::Gamma< 3 > gamma3;
        (*this).rightMultiply(gamma3);
        break;
      }
      case Base::op_GAMMA_4_GAMMA_5:
      {
        Dirac::Gamma< 45 > gamma4gamma5;
        (*this).rightMultiply(gamma4gamma5);
        break;
      }
      default:
        std::cerr << "Error in void Propagator::multiplyOperator(Base::Operator const& O)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        exit(1);
    }
  }


}
