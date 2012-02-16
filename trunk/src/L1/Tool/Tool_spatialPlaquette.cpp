#include "Tool.ih"

double Tool::spatialPlaquette(Core::Field< QCD::Gauge > const &field)
{
  double res = realtr(Path::square(field, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP));
  res += realtr(Path::square(field, Base::idx_X, Base::dir_UP, Base::idx_Z, Base::dir_UP));
  res += realtr(Path::square(field, Base::idx_Y, Base::dir_UP, Base::idx_Z, Base::dir_UP));
  // We define the plaquette as the real part, averaged over the number of distinct orientations and colors
  // So divide by 3 (colors) * 3 (XY, XZ, YZ) = 9
  return (res / 9.0);
}

Core::Correlator< std::complex<double> > Tool::Plaquette_timeslice(Core::Field< QCD::Gauge > const &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo)
{
	Core::Correlator< std::complex<double> > res(  localTrace(Path::square(field, via, dirVia, to, dirTo)));
	
  return (res);
}
