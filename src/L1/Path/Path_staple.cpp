#include "Path.ih"

Core::Field< SU3::Matrix > Path::staple(Core::Field< QCD::Gauge > const &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo)
{
  if (via == to)
    return (dirVia == Base::dir_UP ? const_cast< Core::Field< QCD::Gauge > & >(field).component< SU3::Matrix >(via)
                                   : const_cast< Core::Field< QCD::Gauge > & >(field).component< SU3::Matrix >(via).dagger());
  Core::Field< SU3::Matrix > result(Path::step(field, via, dirVia)); //step1
  Path::step(&result, field, to, dirTo); //step2
  Path::step(&result, field, via, Base::opposite(dirVia)); //step3 (opposite direction)

  result.shift(to, Base::opposite(dirTo)); //back to starting point
  return result;
}


Core::Field< SU3::Matrix > Path::HYPstaple(Core::Field< QCD::Gauge > const &fieldVia,
                                           Core::Field< QCD::Gauge > const &fieldTo,
                                           Base::SpaceTimeIndex via, Base::Direction dirVia,
                                           Base::SpaceTimeIndex to, Base::Direction dirTo)
{
  if (via == to)
    return (dirVia == Base::dir_UP ? const_cast< Core::Field< QCD::Gauge > & >(fieldVia).component< SU3::Matrix >(via)
                                   : const_cast< Core::Field< QCD::Gauge > & >(fieldVia).component< SU3::Matrix >(via).dagger());
  Core::Field< SU3::Matrix > result(Path::step(fieldVia, via, dirVia)); //step1
  Path::step(&result, fieldTo, to, dirTo); //step2
  Path::step(&result, fieldVia, via, Base::opposite(dirVia)); //step3 (opposite direction)

  result.shift(to, Base::opposite(dirTo)); //back to starting point
  return result;
}
