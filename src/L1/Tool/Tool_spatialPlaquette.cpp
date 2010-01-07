#include "Tool.ih"

double Tool::spatialPlaquette(Core::Field< QCD::Gauge > &field)
{
  double res = realtr(Path::square(field, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP));
  res += realtr(Path::square(field, Base::idx_X, Base::dir_UP, Base::idx_Z, Base::dir_UP));
  res += realtr(Path::square(field, Base::idx_Y, Base::dir_UP, Base::idx_Z, Base::dir_UP));
  // We define the plaquette as the real part, averaged over the number of distinct orientations and colors
  // So divide by 3 (colors) * 3 (XY, XZ, YZ) = 9
  return (res / 9.0);
}
