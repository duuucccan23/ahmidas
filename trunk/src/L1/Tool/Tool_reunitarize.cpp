#include "Tool.ih"

void Tool::reunitarize(Core::Field< QCD::Gauge > *field)
{
  field->isolate();
  for (Core::Field< QCD::Gauge >::iterator iter = field->begin(); iter != field->end(); ++iter)
    iter->reunitarize();
}
