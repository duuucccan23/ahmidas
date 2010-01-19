#include "Meson.ih"

namespace Contract
{

//   MesonContraction::MesonContraction()
//   {}
// 
//   MesonContraction::~MesonContraction()
//   {}


  template< size_t IndexSrc, size_t IndexSnk >
  Core::Correlator const &light_meson_twopoint(Core::Propagator const &u, Core::Propagator const &d,
                                                                 Dirac::Gamma< IndexSrc > const &interpolSrc,
                                                                 Dirac::Gamma< IndexSnk > const &interpolSnk)

  {
    assert(u.T() == d.T());
    assert(u.L() == d.L());

    Core::Correlator twopoint(u.L(), u.T());

    std::cerr << "This is about to be implemented!" << std::endl;
    exit(1);

    return twopoint;
  }


  template< size_t IndexSrc, size_t IndexSnk >
  Core::Correlator const &light_meson_twopoint(Core::Propagator const &u, Core::Propagator const &d,
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
