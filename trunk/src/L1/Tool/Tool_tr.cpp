#include "Tool.ih"

std::complex< double > Tool::tr(Core::Field < SU3::Matrix > const &field)
{
  std::complex< double > result(0.0, 0.0);
  for (Core::Field< SU3::Matrix >::const_iterator iter = field.begin(); iter != field.end(); ++iter)
    result += iter->tr();
  size_t vol = field.volume();
  std::complex < double > temp = field.weave(Base::wea_SUM, result);
  temp /= vol;
  return temp;
}
