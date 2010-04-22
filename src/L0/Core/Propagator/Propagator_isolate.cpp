#include "Propagator.ih"

namespace Core
{
  void Propagator::isolate()
  {
    if (*d_references == 1)
      return;

    assert(*d_references > 1);

    *d_references -= 1;
    d_references = new size_t(1);

    Field< QCD::Tensor > *d_components_new = new Field< QCD::Tensor > (L(), T());
    Field< QCD::Tensor >::const_iterator itOld = d_components->begin();
    Field< QCD::Tensor >::iterator itNew = d_components_new->begin();
    while(itOld != d_components->end())
    {
      *itNew = *itOld;
      ++itOld;
      ++itNew;
    }
    d_components = d_components_new;

  }
}
