#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>


namespace Core
{
  // applies the (twisted Mass) Dirac operator with uniform boundary conditions and also its hermitian  conjugate,
  // in order to reconstruct the propagators of the doublet from a combined propagator of the kind (DD^dagger)^-1

  void Propagator::reconstruct_doublet(Propagator &propPlus, Propagator &propMinus,
                                         Field< QCD::Gauge > const &gauge_field,
                                         double const kappa, double const mu, double const thetaT,
                                         double const thetaX, double const thetaY, double const thetaZ) const
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

    // temporal links pick up a factor exp(i*theta*pi/T)
    // spatial links pick up a factor exp(i*theta*pi/L)
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


    propPlus *= std::complex<double>(0.0, 0.0);

    {
      Dirac::Gamma< 4 > gamma4;
      // negative t-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma4*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_T).dagger());
      neighbors.shift(Base::idx_T, Base::dir_UP);
      propPlus += neighbors;
      // poitive t-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma4*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_T, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_T));
      propPlus += neighbors;
    }

    {
      Dirac::Gamma< 1 > gamma1;
      // negative x-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma1*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_X).dagger());
      neighbors.shift(Base::idx_X, Base::dir_UP);
      propPlus += neighbors;
      // positive x-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma1*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_X, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_X));
      propPlus += neighbors;
    }

    {
      Dirac::Gamma< 2 > gamma2;
      // negative y-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma2*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Y).dagger());
      neighbors.shift(Base::idx_Y, Base::dir_UP);
      propPlus += neighbors;
      // positive y-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma2*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_Y, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Y));
      propPlus += neighbors;
    }

   {
      Dirac::Gamma< 3 > gamma3;
      // negative z-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma3*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Z).dagger());
      neighbors.shift(Base::idx_Z, Base::dir_UP);
      propPlus += neighbors;
      // positive z-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma3*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_Z, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply(gauge_field_copy.component< SU3::Matrix >(Base::idx_Z));
      propPlus += neighbors;
    }

    propPlus *= std::complex<double>(-0.5, 0.0);


    // diagonal contribution
    Dirac::Gamma< 5 > gamma5;
    {
      Propagator tmp(*this);
      tmp *= std::complex<double>(0.5/kappa, 0.0);
      propPlus += tmp;
      // twisted mass term ~ i*mu*gamma5
      tmp = gamma5*(*this);
      tmp *= std::complex<double>(0.0, mu);
      propMinus = propPlus; // at this point the calculations for the "plus" and "minus" propagators differ
      propPlus += tmp;
      propMinus -= tmp;
    }

    propPlus.rightMultiply(gamma5);
    propMinus.rightMultiply(gamma5);
    weave.barrier();
  }
}

