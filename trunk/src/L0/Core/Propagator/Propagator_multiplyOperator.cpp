#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

#define __Hermitian_basis__

namespace Core
{


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

#ifdef __Hermitian_basis__

  void Propagator::rightMultiplyOperator(Base::HermitianBilinearOperator const O, bool const twisted_basis)
  {

    Dirac::Gamma< 1 > gamma1;
    Dirac::Gamma< 2 > gamma2;
    Dirac::Gamma< 3 > gamma3;
    Dirac::Gamma< 4 > gamma4;
    Dirac::Gamma< 5 > gamma5;
    std::complex<double > const minus_i(0.0, -1.0);
    std::complex<double > const plus_i(0.0, +1.0);
    isolate();
    switch (O)
    {
      case Base::op_G_0:
        if (twisted_basis)
        {
          // identity (up to a factor)
          (*this)*= plus_i;
        }
        else
        {
          (*this).rightMultiply(gamma5);
        }
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
        (*this).rightMultiply(gamma5);
        (*this).rightMultiply(gamma4);
        (*this)*= minus_i;
        break;
      case Base::op_G_5:
        if (twisted_basis)
        {
          (*this).rightMultiply(gamma1);
          (*this).rightMultiply(gamma4);
          (*this).rightMultiply(gamma5);
        }
        else
        {
          (*this).rightMultiply(gamma1);
          (*this).rightMultiply(gamma4);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_6:
        if (twisted_basis)
        {
          (*this).rightMultiply(gamma2);
          (*this).rightMultiply(gamma4);
          (*this).rightMultiply(gamma5);
        }
        else
        {
          (*this).rightMultiply(gamma2);
          (*this).rightMultiply(gamma4);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_7:
        if (twisted_basis)
        {
          (*this).rightMultiply(gamma3);
          (*this).rightMultiply(gamma4);
          (*this).rightMultiply(gamma5);
        }
        else
        {
          (*this).rightMultiply(gamma3);
          (*this).rightMultiply(gamma4);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_8:
        if (twisted_basis)
        {
          (*this).rightMultiply(gamma5);
          (*this)*= plus_i;
        }
        else
        {
          //unity
        }
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
        if (twisted_basis)
        {
          (*this).rightMultiply(gamma1);
          (*this).rightMultiply(gamma4);
        }
        else
        {
          (*this).rightMultiply(gamma1);
          (*this).rightMultiply(gamma4);
          (*this).rightMultiply(gamma5);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_14:
        if (twisted_basis)
        {
          (*this).rightMultiply(gamma2);
          (*this).rightMultiply(gamma4);
        }
        else
        {
          (*this).rightMultiply(gamma2);
          (*this).rightMultiply(gamma4);
          (*this).rightMultiply(gamma5);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_15:
        if (twisted_basis)
        {
          (*this).rightMultiply(gamma3);
          (*this).rightMultiply(gamma4);
        }
        else
        {
          (*this).rightMultiply(gamma3);
          (*this).rightMultiply(gamma4);
          (*this).rightMultiply(gamma5);
          (*this)*= minus_i;
        }
        break;
      default:
        std::cerr << "Error in void Propagator::multiplyOperator(Base::Operator const& O)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        exit(1);
    }
  }



  void Propagator::leftMultiplyOperator(Base::HermitianBilinearOperator const O, bool const twisted_basis)
  {

    Dirac::Gamma< 1 > gamma1;
    Dirac::Gamma< 2 > gamma2;
    Dirac::Gamma< 3 > gamma3;
    Dirac::Gamma< 4 > gamma4;
    Dirac::Gamma< 5 > gamma5;
    std::complex<double > const minus_i(0.0, -1.0);
    std::complex<double > const plus_i(0.0, +1.0);
    isolate();
    switch (O)
    {
      case Base::op_G_0:
        if (twisted_basis)
        {
          // identity (up to a factor)
          (*this)*= plus_i;
        }
        else
        {
          (*this)*=(gamma5);
        }
        break;
      case Base::op_G_1:
        (*this)*=(gamma1);
        break;
      case Base::op_G_2:
        (*this)*=(gamma2);
        break;
      case Base::op_G_3:
        (*this)*=(gamma3);
        break;
      case Base::op_G_4:
        (*this)*=(gamma5);
        (*this)*=(gamma4);
        (*this)*= minus_i;
        break;
      case Base::op_G_5:
        if (twisted_basis)
        {
          (*this)*=(gamma1);
          (*this)*=(gamma4);
          (*this)*=(gamma5);
        }
        else
        {
          (*this)*=(gamma1);
          (*this)*=(gamma4);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_6:
        if (twisted_basis)
        {
          (*this)*=(gamma2);
          (*this)*=(gamma4);
          (*this)*=(gamma5);
        }
        else
        {
          (*this)*=(gamma2);
          (*this)*=(gamma4);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_7:
        if (twisted_basis)
        {
          (*this)*=(gamma3);
          (*this)*=(gamma4);
          (*this)*=(gamma5);
        }
        else
        {
          (*this)*=(gamma3);
          (*this)*=(gamma4);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_8:
        if (twisted_basis)
        {
          (*this)*=(gamma5);
          (*this)*= plus_i;
        }
        else
        {
          //unity
        }
        break;
      case Base::op_G_9:
        (*this)*=(gamma1);
        (*this)*=(gamma5);
        (*this)*= minus_i;
        break;
      case Base::op_G_10:
        (*this)*=(gamma2);
        (*this)*=(gamma5);
        (*this)*= minus_i;
        break;
      case Base::op_G_11:
        (*this)*=(gamma3);
        (*this)*=(gamma5);
        (*this)*= minus_i;
        break;
      case Base::op_G_12:
        (*this)*=(gamma4);
        break;
      case Base::op_G_13:
        if (twisted_basis)
        {
          (*this)*=(gamma1);
          (*this)*=(gamma4);
        }
        else
        {
          (*this)*=(gamma1);
          (*this)*=(gamma4);
          (*this)*=(gamma5);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_14:
        if (twisted_basis)
        {
          (*this)*=(gamma2);
          (*this)*=(gamma4);
        }
        else
        {
          (*this)*=(gamma2);
          (*this)*=(gamma4);
          (*this)*=(gamma5);
          (*this)*= minus_i;
        }
        break;
      case Base::op_G_15:
        if (twisted_basis)
        {
          (*this)*=(gamma3);
          (*this)*=(gamma4);
        }
        else
        {
          (*this)*=(gamma3);
          (*this)*=(gamma4);
          (*this)*=(gamma5);
          (*this)*= minus_i;
        }
        break;
      default:
        std::cerr << "Error in void Propagator::multiplyOperator(Base::Operator const& O)\n";
        std::cerr << "Operator with no. " << O << " not implemented!" << std::endl;
        exit(1);
    }
  }


#endif

#ifdef __Hermitian_basis_disc_cc__

  void Propagator::rightMultiplyOperator(Base::HermitianBilinearOperator const O)
  {

    Dirac::Gamma< 1 > gamma1;
    Dirac::Gamma< 2 > gamma2;
    Dirac::Gamma< 3 > gamma3;
    Dirac::Gamma< 4 > gamma4;
    Dirac::Gamma< 5 > gamma5;
    std::complex<double >minus_i(0,-1);
    std::complex<double >minus_one(-1,0);
    std::complex<double >plus_i(0,1);
    isolate();
    switch (O)
    {
      case Base::op_G_0:
        (*this).rightMultiply(gamma5);
        break;
      case Base::op_G_1:
        (*this).rightMultiply(gamma1);
        (*this) *=minus_one;
        break;
      case Base::op_G_2:
        (*this).rightMultiply(gamma2);
        (*this) *=minus_one;
        break;
      case Base::op_G_3:
        (*this).rightMultiply(gamma3);
        (*this) *=minus_one;
        break;
      case Base::op_G_4:
        (*this).rightMultiply(gamma5);
        (*this).rightMultiply(gamma4);
        (*this)*= minus_i;
        break;
      case Base::op_G_5:
        (*this).rightMultiply(gamma1);
        (*this).rightMultiply(gamma4);
        (*this)*= plus_i;
        break;
      case Base::op_G_6:
        (*this).rightMultiply(gamma2);
        (*this).rightMultiply(gamma4);
        (*this)*= plus_i;
        break;
      case Base::op_G_7:
        (*this).rightMultiply(gamma3);
        (*this).rightMultiply(gamma4);
        (*this)*= plus_i;
        break;
      case Base::op_G_8:
        (*this) *=minus_one;
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
        (*this) *=minus_one;
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

#endif

}
