#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

namespace Core
{
  // applies the (twisted Mass) Dirac operator with uniform boundary conditions

  Propagator Propagator::applyDiracOperator(Field< QCD::Gauge > const &gauge_field, double const kappa, double const m0,
                                            double const mu_sigma, double const mu_delta, Propagator const &secondProp, bool sign) const
  {
    Base::Weave weave(L(), T());
    
    // mass of doublet is splitted via mu_delta
    double const mass = sign ? m0 + mu_delta : m0 - mu_delta; 
    
    std::complex<double> const phaseFactor(exp(std::complex< double > (0.0, M_PI / double(T()))));
    std::complex<double> const phaseFactorReverse(exp(std::complex< double > (0.0, -M_PI / double(T()))));

    size_t localIndex;

    for(Field< QCD::Gauge >::iterator I(const_cast< Field< QCD::Gauge > & >(gauge_field).begin()); 
        I != (const_cast< Field< QCD::Gauge > & >(gauge_field)).end(); ++I)
    {
      (*I) *= phaseFactor;
    }


    Propagator result(L(), T());
    result *= std::complex<double>(0.0, 0.0);

    {
      Dirac::Gamma< 4 > gamma4;
      // negative t-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma4*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_T).dagger());
      neighbors.shift(Base::idx_T, Base::dir_UP);
      result += neighbors;
      // positive t-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma4*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_T, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_T));
      result += neighbors;
    }

    {
      Dirac::Gamma< 1 > gamma1;
      // negative x-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma1*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_X).dagger());
      neighbors.shift(Base::idx_X, Base::dir_UP);
      result += neighbors;
      // positive x-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma1*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_X, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_X));
      result += neighbors;
    }

    {
      Dirac::Gamma< 2 > gamma2;
      // negative y-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma2*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_Y).dagger());
      neighbors.shift(Base::idx_Y, Base::dir_UP);
      result += neighbors;
      // positive y-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma2*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_Y, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_Y));
      result += neighbors;
    }

   {
      Dirac::Gamma< 3 > gamma3;
      // negative z-direction
      Propagator neighbors(*this);
      neighbors.isolate();
      neighbors += gamma3*neighbors;
      // same procedure as in Transport::step(..., Base::dir_UP)
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_Z).dagger());
      neighbors.shift(Base::idx_Z, Base::dir_UP);
      result += neighbors;
      // positive z-direction
      neighbors = (*this);
      neighbors.isolate();
      neighbors -= gamma3*neighbors;
      // same procedure as in Transport::step(..., Base::dir_DOWN)
      neighbors.shift(Base::idx_Z, Base::dir_DOWN);
      (neighbors.d_components)->rightMultiply((const_cast< Field< QCD::Gauge > & >(gauge_field)).component< SU3::Matrix >(Base::idx_Z));
      result += neighbors;
    }

    result *= std::complex<double>(-0.5, 0.0);


    // diagonal contribution
    {
      Dirac::Gamma< 5 > gamma5;
      Propagator tmp(*this);
      tmp *= std::complex<double>(0.5/kappa, 0.0) + std::complex<double>(mass, 0.0);
      result += tmp;
      // twisted mass term ~ i*mu*gamma5*tau_1 (here we need the corresponding Propagator of the doublet)
      tmp = gamma5*(secondProp);
      tmp *= std::complex<double>(0.0, mu_delta);
      result += tmp;
    }

    // we want to leave the gauge field "untouched"
    for(Field< QCD::Gauge >::iterator I(const_cast< Field< QCD::Gauge > & >(gauge_field).begin()); 
        I != (const_cast< Field< QCD::Gauge > & >(gauge_field)).end(); ++I)
    {
      (*I) *= phaseFactorReverse;
    }
    
    weave.barrier();
    return result;
  }
}
