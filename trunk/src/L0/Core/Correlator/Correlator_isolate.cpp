#include "Correlator.ih"

void Core::Correlator::isolate()
{
  if (*d_references == 1)
    return;

  assert(*d_references > 1);

  *d_references -= 1;
  d_references   = new size_t(1);

  d_weave = new Base::Weave(L(), T());

  Dirac::Matrix *sumTimeslice = new Dirac::Matrix(T());
  std::copy(d_sumTimeslice, d_sumTimeslice + T(), sumTimeslice);
  d_sumTimeslice = sumTimeslice;

  Dirac::Matrix *sumTimeslice_global = new Dirac::Matrix(T());
  std::copy(d_sumTimeslice_global, d_sumTimeslice_global + T(), sumTimeslice_global);
  d_sumTimeslice_global = sumTimeslice_global;

  if (d_data == NULL)
    return;

  Core::Field< Dirac::Matrix > *data = new Core::Field < Dirac::Matrix >(L(), T());
  Core::Field< Dirac::Matrix >::iterator itOld = d_data->begin();
  Core::Field< Dirac::Matrix >::iterator itNew = data->begin();

  while(itOld != d_data->end())
  {
    *itNew = *itOld;
    ++itOld;
    ++itNew;
  }
  d_data = data;
}
