#include "Fuzz.ih"

void Smearing::Fuzz:smear(Fields::GaugeField &field) const
{
  if (d_length <= 1) // Or we'll be making 3 pointless copies...
    return;
  
  accumDirection(field, idx_X);
  accumDirection(field, idx_Y);
  accumDirection(field, idx_Z);
}

