#include "Tool.ih"

Core::Field< std::complex< double > > Tool::localTrace(Core::Field< SU3::Matrix > const &field)
{
  Core::Field< std::complex< double > > result(field.L(), field.T());
  size_t idx = 0;
  for (Core::Field< SU3::Matrix >::const_iterator iter = field.begin(); iter != field.end(); ++iter, ++idx)
    result.fastPhysicalIndex(idx) = iter->tr();
  return result;
}
