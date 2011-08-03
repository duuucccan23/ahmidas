#include "Propagator.ih"

namespace Core
{
  void Propagator::isolate()
  {

    // std::cout << "isolate:\nref = " << d_references << ": refCount (Propagator) = " << *d_references << std::endl;
    if (*d_references == 1)
      return;

    assert(*d_references > 1);

    *d_references -= 1;
    d_references = new size_t(1);


    d_components->destroy();
    d_components = new Field< QCD::Tensor > (*d_components);
    (*d_components).isolate();

  }
}
