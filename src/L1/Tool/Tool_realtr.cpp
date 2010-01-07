#include "Tool.ih"

double Tool::realtr(Core::Field< SU3::Matrix > const &field)
{
  double result(0.0);
  for (Core::Field< SU3::Matrix >::const_iterator iter = field.begin(); iter != field.end(); ++iter)
    result += iter->realtr();

  return field.weave(Base::wea_SUM, result) / field.volume();
}
