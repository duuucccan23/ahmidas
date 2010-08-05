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

      for(size_t t = 0; t < T(); t++)
      {
        Dirac::Matrix second = gamma_5 * d_sumTimeslice_global[t];
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0 * d_sumTimeslice_global[t];
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    case Base::proj_1_PLUS_TM:
    {
      for(size_t t = 0; t < T(); t++)
      {
        Dirac::Gamma< 1 > gamma_1;
        Dirac::Gamma< 51 > gamma5_gamma1;
        Dirac::Matrix second = gamma5_gamma1 * d_sumTimeslice_global[t];
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0 * (gamma_1 * d_sumTimeslice_global[t]);
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    }
    case Base::proj_1_MINUS_TM:
    {
      for(size_t t = 0; t < T(); t++)
      {
        Dirac::Gamma< 1 > gamma_1;
        Dirac::Gamma< 51 > gamma5_gamma1;
        Dirac::Matrix second = gamma5_gamma1 * d_sumTimeslice_global[t];
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0 * (gamma_1 * d_sumTimeslice_global[t]);
        d_sumTimeslice_global[t] *= std::complex< double >(-1, 0);
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    }
    case Base::proj_2_PLUS_TM:
    {
      for(size_t t = 0; t < T(); t++)
      {
        Dirac::Gamma< 2 > gamma_2;
        Dirac::Gamma< 52 > gamma5_gamma2;
        Dirac::Matrix second = gamma5_gamma2 * d_sumTimeslice_global[t];
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0 * (gamma_2 * d_sumTimeslice_global[t]);
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    }
    case Base::proj_2_MINUS_TM:
    {
      for(size_t t = 0; t < T(); t++)
      {
        Dirac::Gamma< 2 > gamma_2;
        Dirac::Gamma< 52 > gamma5_gamma2;
        Dirac::Matrix second = gamma5_gamma2 * d_sumTimeslice_global[t];
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0 * (gamma_2 * d_sumTimeslice_global[t]);
        d_sumTimeslice_global[t] *= std::complex< double >(-1, 0);
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    }
    case Base::proj_3_PLUS_TM:
    {
      for(size_t t = 0; t < T(); t++)
      {
        Dirac::Gamma< 3 > gamma_3;
        Dirac::Gamma< 53 > gamma5_gamma3;
        Dirac::Matrix second = gamma5_gamma3 * d_sumTimeslice_global[t];
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0 * (gamma_3 * d_sumTimeslice_global[t]);
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    }
    case Base::proj_3_MINUS_TM:
    {
      for(size_t t = 0; t < T(); t++)
      {
        Dirac::Gamma< 3 > gamma_3;
        Dirac::Gamma< 53 > gamma5_gamma3;
        Dirac::Matrix second = gamma5_gamma3 * d_sumTimeslice_global[t];
        second *= std::complex< double >(0, 1);
        d_sumTimeslice_global[t] = gamma_0 * (gamma_3 * d_sumTimeslice_global[t]);
        d_sumTimeslice_global[t] *= std::complex< double >(-1, 0);
        d_sumTimeslice_global[t] += second;
        d_sumTimeslice_global[t] *= 0.5;
      }
      break;
    }
    case Base::proj_NO_PROJECTOR:
      // nothing to do
      break;
    default:
      std::cerr << "unknown projector in Correlator::operator*=(BaryonPropagatorProjector const &projector)!" << std::endl;
      std::cerr << "Aborting..." << std::endl;
      exit(1);
  }
}
