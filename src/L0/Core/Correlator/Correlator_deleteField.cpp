#include "Correlator.ih"

void Core::Correlator::deleteField()
{
  isolate();
  d_weave->barrier();
  if (d_data != NULL)
  {
    delete d_data;
    d_data = NULL;
  }
}
