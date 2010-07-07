#include "Tool.ih"

double Tool::temporalPlaquette(Core::Field< QCD::Gauge > const &field)
{
  double res = realtr(Path::square(field, Base::idx_X, Base::dir_UP, Base::idx_T, Base::dir_UP));
  res += realtr(Path::square(field, Base::idx_Y, Base::dir_UP, Base::idx_T, Base::dir_UP));
  res += realtr(Path::square(field, Base::idx_Z, Base::dir_UP, Base::idx_T, Base::dir_UP));
  // We define the plaquette as the real part, averaged over the number of distinct orientations and colors
  // So divide by 3 (colors) * 3 (TX, TY, TZ) = 9
  return (res / 9.0);
}
