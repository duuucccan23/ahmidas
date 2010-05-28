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
          std::complex< double > ZERO(0, 0);


          Propagator copy(*this);
          copy *= gamma4;
//           std::cout << "(*this) *= 0" << std::endl;
          (*this) *= ZERO;

          // std::cout << "Propagator tmp(*this)" << std::endl;
          Propagator tmp(copy);
//           std::cout << "tmp.isolate();" << std::endl;
          tmp.isolate();
//           std::cout << "tmp is isolated;" << std::endl;
          // part 1 : times U_mu and shift up
          (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
          tmp.shift(Base::idx_T, Base::dir_UP);
          (*this) += tmp;

//           std::cout << "part 1 calculated!" << std::endl;

          tmp = copy;
          tmp.isolate();
          // part 2 : shift up and times U_mu^dagger
          tmp.shift(Base::idx_T, Base::dir_UP);
          (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
          tmp *= MINUS_ONE;
          (*this) -= tmp;

//           std::cout << "part 2 calculated!" << std::endl;

          tmp = copy;
          tmp.isolate();
          // part 3 : shift down and times U_mu^dagger
          tmp.shift(Base::idx_T, Base::dir_DOWN);
          (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
          tmp *= MINUS_ONE;
          (*this) -= tmp;

//           std::cout << "part 3 calculated!" << std::endl;

          tmp = copy;
          tmp.isolate();
          // part 4 : times U_mu and shift down
          (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
          tmp.shift(Base::idx_T, Base::dir_DOWN);
          (*this) += tmp;

//           std::cout << "part 4 calculated!" << std::endl;

//           std::cout << "end of scope of tmp, tmp will be deleted now" << std::endl;
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
