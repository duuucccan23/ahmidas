#include "Tool.ih"

double Tool::spatialDownPlaquette(Core::Field< QCD::Gauge > const &field)
{
  // Function to aid in testing the shift routines and the Path::square routine.
  double res = realtr(Path::square(field, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_DOWN));
  res += realtr(Path::square(field, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_DOWN));
  res += realtr(Path::square(field, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_DOWN));
  // We define the plaquette as the real part, averaged over the number of distinct orientations and colors
  // So divide by 3 (colors) * 3 (XY, XZ, YZ) = 9
  return (res / 9.0);
}
