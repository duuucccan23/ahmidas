#include "Fuzz.ih"

void Smearing::Fuzz::accumDirection(Fields::GaugeField &field, SpaceTimeIndex idx) const
{
  Fields::GaugeField shifter(field);
  for (size_t ctr = 1; ctr < d_length; ++ctr)
  {
    shifter.shift(idx, dir_UP);
    field.component(idx).rightMultiply(shifter.component(idx));
  }
}
