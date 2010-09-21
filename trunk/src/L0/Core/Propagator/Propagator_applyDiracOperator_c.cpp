#include "Propagator.ih"

#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>


namespace Core
{
  // applies the (twisted Mass) Dirac operator with fixed boundary conditions
  // (temporal gauge links starting at timeslice t_boundary pick up a factor -1)

  // this is for the light doublet in twisted mass formulation: u quark: mu > 0, d quark: mu < 0
  Propagator Propagator::applyDiracOperator(Field< QCD::Gauge > const &gauge_field, double const kappa, double const mu, size_t const t_boundary) const
  {

    Base::Weave weave(L(), T());
    assert (t_boundary < T());

    size_t localIndex;

    // temporal links at timeslice t_boundary pick up a factor -1
    for(size_t idx_Z = 0; idx_Z < L(); idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L(); idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L(); idx_X++)
        {
          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, t_boundary);
          /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
          if (localIndex == weave.localVolume())
            continue;

          (const_cast< Field< QCD::Gauge > & >(gauge_field))[localIndex][Base::idx_T] *= std::complex<double>(-1.0, 0.0);
        }
      }
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
      tmp *= std::complex<double>(0.5/kappa, 0.0);
      result += tmp;
      // twisted mass term ~ i*mu*gamma5
      tmp = gamma5*(*this);
      tmp *= std::complex<double>(0.0, mu);
      result += tmp;
    }



    // we want to leave the gauge field "untouched"
    for(size_t idx_Z = 0; idx_Z < L(); idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L(); idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L(); idx_X++)
        {
          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, t_boundary);
          /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
          if (localIndex == weave.localVolume())
            continue;

          (const_cast< Field< QCD::Gauge > & >(gauge_field))[localIndex][Base::idx_T] *= std::complex<double>(-1.0, 0.0);
        }
      }
    }


    weave.barrier();
    return result;
  }
}
