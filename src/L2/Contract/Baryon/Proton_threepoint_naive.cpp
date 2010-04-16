#include "Baryon.ih"

namespace Contract
{
  std::vector< Core::Correlator > proton_threepoint_naive(Core:: Propagator const &u1, Core:: Propagator const &u2,
                                                          Core:: Propagator const &d, Core::Field < QCD::Gauge > const * const gauge_field,
                                                          std::vector< Base::Operator > const &ops,
                                                          size_t const t_src, size_t const t_snk)
  {
    assert(u1.L() == d.L() && u1.T() == d.T() && u2.L() == d.L() && u2.T() == d.T());
    if(gauge_field != NULL)
      assert(gauge_field.L() == d.L() && gauge_field.T() == d.T);





    for (size_t opIdx=0; opIdx<my_operators.size(); opIdx++)
    {
      Core::Correlator threepoint_UU(u.L(), u.T(), threepoint[2*opIdx  ]);
      Core::Correlator threepoint_DD(u.L(), u.T(), threepoint[2*opIdx+1]);
    }
  }
}
