#ifndef GUARD_CONTRACT_WICK_H
#define GUARD_CONTRACT_WICK_H

#include <L0/Core/Field.h>
#include <L0/Core/Buffer.h>
#include <L0/QCD/Spinor.h>
#include <L1/QCD/Tensor.h>

namespace Contract
{
  Core::Buffer< std::complex, L, T > wick(Core::Field< QCD::Spinor, L, T >, Core);
}

#endif
