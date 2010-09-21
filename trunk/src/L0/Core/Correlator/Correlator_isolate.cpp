#include "Correlator.ih"

template< typename Datatype >
void Core::Correlator< Datatype >::isolate()
{
  if (*d_references == 1)
    return;

  assert(*d_references > 1);

  *d_references -= 1;
  d_references   = new size_t(1);

  d_weave = new Base::Weave(L(), T());

  Dirac::Matrix *sumTimeslice = new Datatype[T()];
//   for(size_t idx_T = 0; idx_T < T(); idx_T++)
//     sumTimeslice[idx_T] = Dirac::Matrix(d_sumTimeslice[idx_T]);
  std::copy(d_sumTimeslice, d_sumTimeslice + T(), sumTimeslice);
  d_sumTimeslice = sumTimeslice;

  Dirac::Matrix *sumTimeslice_global = new Datatype[T()];
//   for(size_t idx_T = 0; idx_T < T(); idx_T++)
//     sumTimeslice[idx_T] = Dirac::Matrix(d_sumTimeslice[idx_T]);
  std::copy(d_sumTimeslice_global, d_sumTimeslice_global + T(), sumTimeslice_global);
  d_sumTimeslice_global = sumTimeslice_global;

  d_weave->barrier();

  if (d_data == NULL)
    return;

  d_data = new Field< Datatype > (*d_data);
  (*d_data).isolate();

  d_weave->barrier();
}
