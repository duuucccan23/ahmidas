#include "Path.ih"

void Path::step(Core::Field< SU3::Matrix > *path, Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex idx, Base::Direction dir, size_t nsteps)
{
  if (dir == Base::dir_UP)
  {
    for ( ; nsteps > 0; nsteps--)
    {
      path->rightMultiply(field.component< SU3::Matrix >(idx));
      path->shift(idx, Base::dir_UP);
    }
    return;
  }
  for ( ; nsteps > 0; nsteps--)
  {
    path->shift(idx, Base::dir_DOWN);
    path->rightMultiply(field.component< SU3::Matrix >(idx).dagger());
  }
}

//To be used in first steps
Core::Field< SU3::Matrix > Path::step(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex idx, Base::Direction dir, size_t nsteps)
{
  if (dir == Base::dir_DOWN)
  {
    if (nsteps == 1)
      return Core::Field< SU3::Matrix > (field.component< SU3::Matrix >(idx).dagger());
    Core::Field< SU3::Matrix > result(field.component< SU3::Matrix >(idx).dagger());
    step(&result, field, idx, dir, --nsteps);
    return result;
  }
  Core::Field< SU3::Matrix > result(field.component< SU3::Matrix >(idx));
  result.shift(idx, Base::dir_UP);
  if (nsteps == 1)
    return result;
  step(&result, field, idx, dir, --nsteps);
  return result;
}
