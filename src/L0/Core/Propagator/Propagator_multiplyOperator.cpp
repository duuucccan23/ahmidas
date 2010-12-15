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
        std::cout << "calculating (symmetrized) conserved vector current, temporal component ... " << std::endl;
        
        assert(gauge_field != NULL);

        Dirac::Gamma< 4 > gamma4;

        {
        Propagator tmp(fromRight);
        std::cout << "part 1 ... ";
        std::cout.flush();
        tmp.isolate();
        // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->rightMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        Propagator copy(*this);
        // gamma_mu - identity
        copy *= gamma4;
        copy -= (*this);
        result = copy.contract(tmp);
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
        Propagator copy(*this);
        // gamma_mu + identity
        copy *= gamma4;
        copy += (*this);
        (*result) += (*(copy.contract(tmp)));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 3 ... ";
        std::cout.flush();
        // gamma_mu + identity
        tmp *= gamma4;
        tmp += (*this);
        // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
        tmp.shift(Base::idx_T, Base::dir_DOWN);
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T).dagger());
        (*result) += *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }

        {
        Propagator tmp(*this);
        tmp.isolate();
        std::cout << "part 4 ... ";
        std::cout.flush();
        // gamma_mu - identity
        tmp *= gamma4;
        tmp -= (*this);
        // part 4 : multiply this by U_mu from right and and shift it up
        (tmp.d_components)->leftMultiply((*gauge_field).component< SU3::Matrix >(Base::idx_T));
        tmp.shift(Base::idx_T, Base::dir_UP);
        (*result) += *(tmp.contract(fromRight));
        std::cout << "done." << std::endl;
        }


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
      case Base::op_GAMMA_15:
      {
        Dirac::Gamma< 15 > gamma1gamma5;
        (*this) *= gamma1gamma5;
        break;
      }
      case Base::op_GAMMA_25:
      {
        Dirac::Gamma< 25 > gamma2gamma5;
        (*this) *= gamma2gamma5;
        break;
      }
      case Base::op_GAMMA_35:
      {
        Dirac::Gamma< 35 > gamma3gamma5;
        (*this) *= gamma3gamma5;
        break;
      }
      case Base::op_GAMMA_45:
      {
        Dirac::Gamma< 45 > gamma4gamma5;
        (*this) *= gamma4gamma5;
        break;
      }
      case Base::op_GAMMA_12:
      {
        Dirac::Gamma< 12 > gamma1gamma2;
        (*this) *= gamma1gamma2;
        break;
      }
      case Base::op_GAMMA_13:
      {
        Dirac::Gamma< 13 > gamma1gamma3;
        (*this) *= gamma1gamma3;
        break;
      }
      case Base::op_GAMMA_14:
      {
        Dirac::Gamma< 14 > gamma1gamma4;
        (*this) *= gamma1gamma4;
        break;
      }
      case Base::op_GAMMA_23:
      {
        Dirac::Gamma< 23 > gamma2gamma3;
        (*this) *= gamma2gamma3;
        break;
      }
      case Base::op_GAMMA_24:
      {
        Dirac::Gamma< 24 > gamma2gamma4;
        (*this) *= gamma2gamma4;
        break;
      }
      case Base::op_GAMMA_34:
      {
        Dirac::Gamma< 34 > gamma3gamma4;
        (*this) *= gamma3gamma4;
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
      case Base::op_GAMMA_5:
      {
        Dirac::Gamma< 5 > gamma5;
        (*this).rightMultiply(gamma5);
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
      case Base::op_GAMMA_15:
      {
        Dirac::Gamma< 15 > gamma1gamma5;
        (*this).rightMultiply(gamma1gamma5);
        break;
      }
      case Base::op_GAMMA_25:
      {
        Dirac::Gamma< 25 > gamma2gamma5;
        (*this).rightMultiply(gamma2gamma5);
        break;
      }
	  case Base::op_GAMMA_35:
	  {
		  Dirac::Gamma< 35 > gamma3gamma5;
		  (*this).rightMultiply(gamma3gamma5);
		  break;
	  }
	  case Base::op_GAMMA_45:
	  {
		  Dirac::Gamma< 45 > gamma4gamma5;
		  (*this).rightMultiply(gamma4gamma5);
		  break;
	  }
	  case Base::op_GAMMA_12:
	  {
		  Dirac::Gamma< 12 > gamma1gamma2;
		  (*this).rightMultiply(gamma1gamma2);
		  break;
	  }
	  case Base::op_GAMMA_13:
	  {
		  Dirac::Gamma< 13 > gamma1gamma3;
		  (*this).rightMultiply(gamma1gamma3);
		  break;
	  }
	  case Base::op_GAMMA_14:
	  {
		  Dirac::Gamma< 14 > gamma1gamma4;
		  (*this).rightMultiply(gamma1gamma4);
		  break;
	  }
	  case Base::op_GAMMA_23:
	  {
		  Dirac::Gamma< 23 > gamma2gamma3;
		  (*this).rightMultiply(gamma2gamma3);
		  break;
	  }
	  case Base::op_GAMMA_24:
	  {
		  Dirac::Gamma< 24 > gamma2gamma4;
		  (*this).rightMultiply(gamma2gamma4);
		  break;
	  }
	  case Base::op_GAMMA_34:
	  {
		  Dirac::Gamma< 34 > gamma3gamma4;
		  (*this).rightMultiply(gamma3gamma4);
		  break;
	  }
	  default:
	  std::cerr << "Error in void Propagator::multiplyOperator(Base::Operator const& O)\n";
	  std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
	  exit(1);
	}
  }

  void Propagator::rightMultiplyOperator(Base::HermitianBilinearOperator const O)
  {

	  Dirac::Gamma< 1 > gamma1;
	  Dirac::Gamma< 2 > gamma2;
	  Dirac::Gamma< 3 > gamma3;
	  Dirac::Gamma< 4 > gamma4;
	  Dirac::Gamma< 5 > gamma5;
	  std::complex<double >minus_i(0,-1); 
	  isolate();
	  switch (O)
	  {
		  case Base::op_G_0:
			  (*this).rightMultiply(gamma5);
			  break;
		  case Base::op_G_1:
			  (*this).rightMultiply(gamma1);
			  break;
		  case Base::op_G_2:
			  (*this).rightMultiply(gamma2);
			  break;
		  case Base::op_G_3:
			  (*this).rightMultiply(gamma3);
			  break;
		  case Base::op_G_4:
			  (*this).rightMultiply(gamma4);
			  (*this).rightMultiply(gamma5);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_5:
			  (*this).rightMultiply(gamma1);
			  (*this).rightMultiply(gamma4);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_6:
			  (*this).rightMultiply(gamma2);
			  (*this).rightMultiply(gamma4);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_7:
			  (*this).rightMultiply(gamma3);
			  (*this).rightMultiply(gamma4);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_8:
			  break;
		  case Base::op_G_9:
			  (*this).rightMultiply(gamma1);
			  (*this).rightMultiply(gamma5);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_10:
			  (*this).rightMultiply(gamma2);
			  (*this).rightMultiply(gamma5);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_11:
			  (*this).rightMultiply(gamma3);
			  (*this).rightMultiply(gamma5);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_12:
			  (*this).rightMultiply(gamma4);
			  break;
		  case Base::op_G_13:
			  (*this).rightMultiply(gamma1);
			  (*this).rightMultiply(gamma4);
			  (*this).rightMultiply(gamma5);
			  (*this)*= minus_i;
			  break;
		  case Base::op_G_14:
			  (*this).rightMultiply(gamma2);
			  (*this).rightMultiply(gamma4);
			  (*this).rightMultiply(gamma5);
			  (*this)*= minus_i;

			  break;
		  case Base::op_G_15:
			  (*this).rightMultiply(gamma3);
			  (*this).rightMultiply(gamma4);
			  (*this).rightMultiply(gamma5);
			  (*this)*= minus_i;

			  break;
		  default:
			  std::cerr << "Error in void Propagator::multiplyOperator(Base::Operator const& O)\n";
			  std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
			  exit(1);
	  }
  }


}
