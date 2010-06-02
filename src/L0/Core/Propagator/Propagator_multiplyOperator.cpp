#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

namespace Core
{

  void Propagator::multiplyOperator(Base::Operator const O, Field< QCD::Gauge > * const gauge_field)
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
      case Base::op_O44:
      {

        std::cout << "calculating product with operator O44 ... " << std::endl;
        assert(gauge_field != NULL);

        Dirac::Gamma< 4 > gamma4;
        std::complex< double > ZERO(0, 0);

        Propagator copy(*this);
        copy *= gamma4;
        (*this) *= ZERO;

        {
        Propagator tmp(copy);
        std::cout << "part 1 ... ";
        std::cout.flush();
        tmp.isolate();
        // part 1 : times U_mu and shift up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        (*this) += tmp;
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(copy);
        tmp.isolate();
        std::cout << "part 2 ... ";
        std::cout.flush();
        // part 2 : shift up and times U_mu^dagger
        tmp.shift(Base::idx_T, Base::dir_UP);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        (*this) -= tmp;
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(copy);
        tmp.isolate();
        std::cout << "part 3 ... ";
        std::cout.flush();
        // part 3 : shift down and times U_mu^dagger
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        (*this) -= tmp;
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(copy);
        tmp.isolate();
        std::cout << "part 4 ... ";
        std::cout.flush();
        // part 4 : times U_mu and shift down
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (*this) += tmp;
        std::cout << "done." << std::endl;
        }
        break;
      }
      case Base::op_O11:
      {
        assert(gauge_field != NULL);

        Dirac::Gamma< 1 > gamma1;
        std::complex< double > ZERO(0, 0);

        Propagator copy(*this);
        copy *= gamma1;
        (*this) *= ZERO;
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 1 : times U_mu and shift up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_UP);
        (*this) += tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 2 : shift up and times U_mu^dagger
        tmp.shift(Base::idx_X, Base::dir_UP);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        (*this) -= tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 3 : shift down and times U_mu^dagger
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X).dagger());
        (*this) -= tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 4 : times U_mu and shift down
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_X));
        tmp.shift(Base::idx_X, Base::dir_DOWN);
        (*this) += tmp;
        }

        break;
      }
      case Base::op_O22:
      {
        assert(gauge_field != NULL);

        Dirac::Gamma< 2 > gamma2;
        std::complex< double > ZERO(0, 0);

        Propagator copy(*this);
        copy *= gamma2;
        (*this) *= ZERO;

        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 1 : times U_mu and shift up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_UP);
        (*this) += tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 2 : shift up and times U_mu^dagger
        tmp.shift(Base::idx_Y, Base::dir_UP);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        (*this) -= tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 3 : shift down and times U_mu^dagger
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y).dagger());
        (*this) -= tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 4 : times U_mu and shift down
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Y));
        tmp.shift(Base::idx_Y, Base::dir_DOWN);
        (*this) += tmp;
        }
        break;
      }
      case Base::op_O33:
      {
        assert(gauge_field != NULL);

        Dirac::Gamma< 3 > gamma3;
        std::complex< double > ZERO(0, 0);

        Propagator copy(*this);
        copy *= gamma3;
        (*this) *= ZERO;
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 1 : times U_mu and shift up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_UP);
        (*this) += tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 2 : shift up and times U_mu^dagger
        tmp.shift(Base::idx_Z, Base::dir_UP);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        (*this) -= tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 3 : shift down and times U_mu^dagger
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z).dagger());
        (*this) -= tmp;
        }
        {
        Propagator tmp(copy);
        tmp.isolate();
        // part 4 : times U_mu and shift down
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_Z));
        tmp.shift(Base::idx_Z, Base::dir_DOWN);
        (*this) += tmp;
        }

        break;
      }
      default:
        std::cerr << "Error in void Propagator::operator*=(Base::Operator const& O)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        exit(1);
    }
  }

}
