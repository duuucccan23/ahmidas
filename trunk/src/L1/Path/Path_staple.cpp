#include "Path.ih"

Core::Field< SU3::Matrix > Path::staple(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo)
{
  if (via == to)
    return (dirVia == Base::dir_UP ? field.component< SU3::Matrix >(via) : field.component< SU3::Matrix >(via).dagger());
  Core::Field< SU3::Matrix > result(Path::step(field, via, dirVia)); //step1
  Path::step(&result, field, to, dirTo); //step2
  Path::step(&result, field, via, opposite(dirVia)); //step3 (opposite direction)
  result.shift(to, opposite(dirTo)); //back to starting point
  return result;
}
