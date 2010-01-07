#include "Tool.ih"

Core::Field< double > Tool::localRealTrace(Core::Field< SU3::Matrix > const &field)
{
  Core::Field< double > result(field.L(), field.T());
  size_t idx = 0;
  for (Core::Field< SU3::Matrix >::const_iterator iter = field.begin(); iter != field.end(); ++iter, ++idx)
    result.fastPhysicalIndex(idx) = iter->realtr();
  return result;
}
