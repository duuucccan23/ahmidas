#include "Path.ih"

void Path::step(Core::Field< SU3::Matrix > *path, Core::Field< QCD::Gauge > const &field, Base::SpaceTimeIndex idx, Base::Direction dir, size_t nsteps)
{
  if (dir == Base::dir_UP)
  {
    for ( ; nsteps > 0; nsteps--)
    {
      path->rightMultiply(const_cast< Core::Field< QCD::Gauge > & >(field).component< SU3::Matrix >(idx));
      // send data from  site x to site x+mu
      path->shift(idx, Base::dir_UP);
    }
    return;
  }
  for ( ; nsteps > 0; nsteps--)
  {
    // send data from  site x to site x-mu
    path->shift(idx, Base::dir_DOWN);
    path->rightMultiply(const_cast< Core::Field< QCD::Gauge > & >(field).component< SU3::Matrix >(idx).dagger());
  }
}

//To be used in first steps
Core::Field< SU3::Matrix > Path::step(Core::Field< QCD::Gauge > const &field, Base::SpaceTimeIndex idx, Base::Direction dir, size_t nsteps)
{
  if (dir == Base::dir_DOWN)
  {
    //field.shift(idx, Base::dir_DOWN);
    Core::Field< SU3::Matrix > result(const_cast< Core::Field< QCD::Gauge > & >(field).component< SU3::Matrix >(idx).dagger());
    //field.shift(idx, Base::dir_UP);
    //result.shift(idx, Base::dir_DOWN);
    if (nsteps == 1)
      return result;
    step(&result, field, idx, dir, --nsteps);
    return result;
  }
  else
  {
    Core::Field< SU3::Matrix > result(const_cast< Core::Field< QCD::Gauge > & >(field).component< SU3::Matrix >(idx));
    result.isolate();
    result.shift(idx, Base::dir_UP);
    if (nsteps == 1)
      return result;
    step(&result, field, idx, dir, --nsteps);
    return result;
  }
}
