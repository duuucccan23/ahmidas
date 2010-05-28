#pragma once

#include <complex>
#include <iostream>
#include <iomanip>

namespace Dirac
{
  class Identity
  {};

  std::ostream &operator<<(std::ostream &out, Dirac::Identity const &id);
}

#include "Identity/Identity_cout_operator_lshift.template"
