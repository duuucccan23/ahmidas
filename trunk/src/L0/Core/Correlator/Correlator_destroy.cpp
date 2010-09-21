#include "Correlator.ih"

template< typename Datatype >
void Core::Correlator< Datatype >::destroy()
{
  assert(*d_references > 0);
  *d_references -= 1;
  if (*d_references == 0)
  {
    delete    d_references;
    delete    d_weave;
    delete [] d_sumTimeslice;
    delete [] d_sumTimeslice_global;
    if (d_data != NULL)
      delete d_data;
  }
  else if (d_data != NULL)
  {
    (*d_data).destroy(); // for proper reference counting of the Field
  }
}
