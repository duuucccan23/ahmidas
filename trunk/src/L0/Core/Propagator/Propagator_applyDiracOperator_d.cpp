#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>


namespace Core
{
  // applies the (twisted Mass) Dirac operator with uniform or fixed boundary conditions
  // This member function can also apply only a part of the Dirac operator (usefull at least when computing disconnected diagrams with the variance reduction based on the Hopping parameter expansion.
  // The function can apply D_OS, A, H or B and Bdagger. The operator are defined as : (D_OS standard OS operator)
  // M = 2  kappa D_OS = ( A + H) )
  // where A is the diagonal part of D_OS and H is the hopping term.
  // B = 1/A
  // It is then straightforward to show that Bdagger = A/(1+4kappa^2 mu^2)
  // this is for the light doublet in twisted mass formulation: u quark: mu > 0, d quark: mu < 0
  Propagator Propagator::applyDiracOperator(Field< QCD::Gauge > const &gauge_field, double const kappa, double const mu,
                                            double const thetaT, double const thetaX, double const thetaY, double const thetaZ,Base::DiracOperator const O) const
  {

	  Base::Weave weave(L(), T());

	  std::complex<double> const phaseFactorT(exp(std::complex< double > (0.0, thetaT * M_PI / double(T()))));
	  std::complex<double> const phaseFactorReverseT(exp(std::complex< double > (0.0, -thetaT * M_PI / double(T()))));
	  std::complex<double> const phaseFactorX(exp(std::complex< double > (0.0, thetaX * M_PI / double(L()))));
	  std::complex<double> const phaseFactorReverseX(exp(std::complex< double > (0.0, -thetaX * M_PI / double(L()))));
	  std::complex<double> const phaseFactorY(exp(std::complex< double > (0.0, thetaY * M_PI / double(L()))));
	  std::complex<double> const phaseFactorReverseY(exp(std::complex< double > (0.0, -thetaY * M_PI / double(L()))));
	  std::complex<double> const phaseFactorZ(exp(std::complex< double > (0.0, thetaZ * M_PI / double(L()))));
	  std::complex<double> const phaseFactorReverseZ(exp(std::complex< double > (0.0, -thetaZ * M_PI / double(L()))));

	  // this is more memory consuming, but safer than modifying the original (const) gauge field with a const_cast
	  Field< QCD::Gauge > gauge_field_copy(gauge_field);
	  gauge_field_copy.isolate();

	  for(Field< QCD::Gauge >::iterator I(gauge_field_copy.begin()); I != gauge_field_copy.end(); ++I)
	  {
		  if(thetaT != 0.0)
			  (*I)[Base::idx_T] *= phaseFactorT;
		  if(thetaX != 0.0)
			  (*I)[Base::idx_X] *= phaseFactorX;
		  if(thetaY != 0.0)
			  (*I)[Base::idx_Y] *= phaseFactorY;
		  if(thetaZ != 0.0)
			  (*I)[Base::idx_Z] *= phaseFactorZ;
	  }

	  //a normalization constant
	  double const nrm=1./(1+4*kappa*kappa*mu*mu);
	  
	  Propagator result(L(), T());
	  result *= std::complex<double>(0.0, 0.0);

	  switch(O)
	  {
		  case Base::Full:
			  {
				  {
					  Dirac::Gamma< 4 > gamma4;
					  // negative t-direction
					  Propagator neighbors(*this);
					  //  neighbors.isolate();
					  neighbors += gamma4*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_T).dagger());
					  neighbors.shift(Base::idx_T, Base::dir_UP);
					  result += neighbors;
					  // positive t-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma4*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_T, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_T));
					  result += neighbors;
				  }

				  {
					  Dirac::Gamma< 1 > gamma1;
					  // negative x-direction
					  Propagator neighbors(*this);
					  //neighbors.isolate();
					  neighbors += gamma1*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_X).dagger());
					  neighbors.shift(Base::idx_X, Base::dir_UP);
					  result += neighbors;
					  // positive x-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma1*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_X, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_X));
					  result += neighbors;
				  }

				  {
					  Dirac::Gamma< 2 > gamma2;
					  // negative y-direction
					  Propagator neighbors(*this);
					  //neighbors.isolate();
					  neighbors += gamma2*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Y).dagger());
					  neighbors.shift(Base::idx_Y, Base::dir_UP);
					  result += neighbors;
					  // positive y-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma2*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_Y, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Y));
					  result += neighbors;
				  }

				  {
					  Dirac::Gamma< 3 > gamma3;
					  // negative z-direction
					  Propagator neighbors(*this);
					  //neighbors.isolate();
					  neighbors += gamma3*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Z).dagger());
					  neighbors.shift(Base::idx_Z, Base::dir_UP);
					  result += neighbors;
					  // positive z-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma3*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_Z, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Z));
					  result += neighbors;
				  }

				  result *= std::complex<double>(-0.5, 0.0);

				  // diagonal contribution
				  {
					  Dirac::Gamma< 5 > gamma5;
					  Propagator tmp(*this);
					  tmp *= std::complex<double>(0.5/kappa, 0.0);
					  result += tmp;
					  // twisted mass term ~ i*mu*gamma5
					  tmp = gamma5*(*this);
					  tmp *= std::complex<double>(0.0, mu);
					  result += tmp;
				  }


				  break;
			  }

		  case Base::H:
			  {

				  {
					  Dirac::Gamma< 4 > gamma4;
					  // negative t-direction
					  Propagator neighbors(*this);
					  //neighbors.isolate();
					  neighbors += gamma4*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_T).dagger());
					  neighbors.shift(Base::idx_T, Base::dir_UP);
					  result += neighbors;
					  // positive t-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma4*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_T, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_T));
					  result += neighbors;
				  }

				  {
					  Dirac::Gamma< 1 > gamma1;
					  // negative x-direction
					  Propagator neighbors(*this);
					  //neighbors.isolate();
					  neighbors += gamma1*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_X).dagger());
					  neighbors.shift(Base::idx_X, Base::dir_UP);
					  result += neighbors;
					  // positive x-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma1*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_X, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_X));
					  result += neighbors;
				  }

				  {
					  Dirac::Gamma< 2 > gamma2;
					  // negative y-direction
					  Propagator neighbors(*this);
					  //neighbors.isolate();
					  neighbors += gamma2*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Y).dagger());
					  neighbors.shift(Base::idx_Y, Base::dir_UP);
					  result += neighbors;
					  // positive y-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma2*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_Y, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Y));
					  result += neighbors;
				  }

				  {
					  Dirac::Gamma< 3 > gamma3;
					  // negative z-direction
					  Propagator neighbors(*this);
					  //neighbors.isolate();
					  neighbors += gamma3*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_UP)
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Z).dagger());
					  neighbors.shift(Base::idx_Z, Base::dir_UP);
					  result += neighbors;
					  // positive z-direction
					  neighbors = (*this);
					  //neighbors.isolate();
					  neighbors -= gamma3*neighbors;
					  // same procedure as in Transport::step(..., Base::dir_DOWN)
					  neighbors.shift(Base::idx_Z, Base::dir_DOWN);
					  (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Z));
					  result += neighbors;
				  }

				  result *= std::complex<double>(-kappa, 0.0);
			
				  break;
			  }



		  case Base::A:
			  {
				  Dirac::Gamma< 5 > gamma5;
				  Propagator tmp(*this);
				  tmp *= std::complex<double>(1, 0.0);
				  result += tmp;
				  // twisted mass term ~ 2*kappa*i*mu*gamma5
				  tmp = gamma5*(*this);
				  tmp *= std::complex<double>(0.0,2*kappa* mu);
				  result += tmp;
				  break;
			  }

		  case Base::B:
			  {

				  Dirac::Gamma< 5 > gamma5;
				  Propagator tmp(*this);
				  tmp *= std::complex<double>(1, 0.0);
				  result += tmp;
				  // twisted mass term ~ 2 kappa * i*mu*gamma5
				  tmp = gamma5*(*this);
				  tmp *= std::complex<double>(0.0,-2*kappa* mu);
				  tmp *= nrm;
				  result += tmp;
				  break;
			  }
		  case Base::Bdagger:
			  {
				  Dirac::Gamma< 5 > gamma5;
				  Propagator tmp(*this);
				  tmp *= std::complex<double>(1, 0.0);
				  result += tmp;
				  // twisted mass term ~ 2 kappa i*mu*gamma5
				  tmp = gamma5*(*this);
				  tmp *= std::complex<double>(0.0,+2*kappa* mu);
				  tmp *= nrm;
				  result += tmp;

				  break;
			  }

	  };


	  weave.barrier();
	  return result;
  }
}
