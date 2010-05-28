#include "Propagator.ih"

namespace Core
{

  void Propagator::destroy()
  {
    // std::cout << "ref = " << d_references << ": refCount (Propagator) = " << *d_references << std::endl;

    assert(*d_references > 0);

    *d_references -= 1;

    if (*d_references > 0)
    {
      (*d_components).destroy(); // for proper reference counting of the Field
    }
    else if (*d_references == 0)
    {
      delete d_components;
      delete d_references;
    }
  }
}
