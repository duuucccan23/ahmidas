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
        {
          Dirac::Gamma< 4 > gamma4;
          (*this) *= gamma4;
        }
        break;
      }
      case Base::op_GAMMA_1:
      {
        {
          Dirac::Gamma< 1 > gamma1;
          (*this) *= gamma1;
        }
        break;
      }
      case Base::op_GAMMA_2:
      {
        {
          Dirac::Gamma< 2 > gamma2;
          (*this) *= gamma2;
        }
        break;
      }
      case Base::op_GAMMA_3:
      {
        {
          Dirac::Gamma< 3 > gamma3;
          (*this) *= gamma3;
        }
        break;
      }
      case Base::op_O44:
      {
        {
          assert(gauge_field != NULL);

          Dirac::Gamma< 4 > gamma4;
          std::complex< double > MINUS_ONE(-1, 0);

          (*this) *= gamma4;

          Propagator result(*this);
          // part 1 : times U_mu and shift up
          result.isolate();
          (result.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
          result.shift(Base::idx_T, Base::dir_UP);

          Propagator tmp(*this);
          tmp.isolate();
          // part 2 : shift up and times U_mu^dagger
          tmp.shift(Base::idx_T, Base::dir_UP);
          (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
          tmp *= MINUS_ONE;
          (*(result.d_components)) += (*(tmp.d_components));

          tmp = (*this);
          // part 3 : shift down and times U_mu^dagger
          tmp.shift(Base::idx_T, Base::dir_DOWN);
          (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
          tmp *= MINUS_ONE;
          (*(result.d_components)) += (*(tmp.d_components));

          tmp = (*this);
          // part 4 : times U_mu and shift down
          (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
          tmp.shift(Base::idx_T, Base::dir_DOWN);
          (*(result.d_components)) += (*(tmp.d_components));
          (*this) = Propagator(result);
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
