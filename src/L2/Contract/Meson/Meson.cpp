#include "Meson.ih"

namespace Contract
{


  template< size_t IndexSrc, size_t IndexSnk >
  Core::Correlator light_meson_twopoint(Core::Propagator const &u, Core::Propagator const &d,
                                                                 Dirac::Gamma< IndexSrc > const &interpolSrc,
                                                                 Dirac::Gamma< IndexSnk > const &interpolSnk,
                                                                 double const *momentum)
  {
    assert(u.T() == d.T());
    assert(u.L() == d.L());
    Core::Correlator twopoint(u.L(), u.T());

    std::cerr << "This has not been implemented yet!" << std::endl;
    exit(1);

    return twopoint;
  }

}
