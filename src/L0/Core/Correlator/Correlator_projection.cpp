#include "Correlator.ih"

// this has only effect on d_sumTimeslice (maybe this has to be changed in the future;
void Core::Correlator::operator*=(Base::BaryonPropagatorProjector const projector)
{
  assert(d_sumTimeslice_global != NULL);
  Dirac::Gamma< 5 > gamma_5;
  Dirac::Gamma< 4 > gamma_0;
  isolate();
  switch (projector)
  {
    case Base::proj_PARITY_PLUS_TM:
      // projector is gamma0 + i*gamma5

      for(size_t t=0; t<T; t++)
      {
        QCD::reducedTensor second = gamma_5*d_sumTimeslice_global[t];
        //second *= gamma_5;
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0*d_sumTimeslice_global[t];
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    default:
      std::cerr << "unknown projector in Correlator::operator*=(BaryonPropagatorProjector const &projector)!" << std::endl;
      std::cerr << "Aborting..." << std::endl;
      exit(1);
  }
}
