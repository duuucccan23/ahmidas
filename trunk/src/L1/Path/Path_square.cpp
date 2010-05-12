#include "Path.ih"

Core::Field< SU3::Matrix > Path::square(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo)
{ //Note that we shift the result field around, it will be smaller
  if (via == to)
    return Core::Field< SU3::Matrix > (SU3::Matrix::identity(), field.L(), field.T());
  Core::Field< SU3::Matrix > result(Path::step(field, via, dirVia)); //step1

  Path::step(&result, field, to, dirTo); //step2
  Path::step(&result, field, via, Base::opposite(dirVia)); //step3
  Path::step(&result, field, to,  Base::opposite(dirTo)); //back to starting point
  return result;
}
