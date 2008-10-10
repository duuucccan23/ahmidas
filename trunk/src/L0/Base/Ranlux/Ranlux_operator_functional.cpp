#include "Ranlux.ih"

double Base::Ranlux::operator()()
{
  for (size_t ctr = 0; ctr < d_lux; ++ctr)
  {
    ++d_index;
    d_index %= d_period;
    d_state[d_index] -= d_state[(d_index + d_delta_s) % d_period] + d_carry;
    d_carry = d_state[d_index] < 0.0 ? d_delta_b : 0.0;
    d_state[d_index] += d_state[d_index] < 0.0 ? 1.0 : 0.0;
  }
  return d_state[d_index];
}
