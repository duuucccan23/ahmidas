#include "Tool.ih"

void Tool::fixCoulombGauge(Core::Field< QCD::Gauge > *field)
{
  // We use the method described in Phys. Rev. D37, p. 1581. (Davies et al.)
  // It exists of a steepest descend, followed by reunitarization.
  // Currently, no Fourier acceleration is being used.
  double const alpha = 0.05;

  Core::Field< SU3::Matrix > transform(field->L(), field->T());
  Core::Field< SU3::Matrix > scratch(field->L(), field->T());
  Core::Field< SU3::Matrix > scratch2(field->L(), field->T());

  std::complex< double > gaugeTerm(1.0, 0.0);
  size_t iterations;

  for (iterations = 0; iterations < 1000; ++iterations)
  {
    setToIdentity(&transform);

    // Dir X
    scratch  = field->component< SU3::Matrix >(Base::idx_X).dagger(); // + U^\dagger_\nu(x)
    scratch -= field->component< SU3::Matrix >(Base::idx_X);          // - U_\nu(x)

    field->shift(Base::idx_X, Base::dir_UP);
    scratch -= field->component< SU3::Matrix >(Base::idx_X).dagger(); // - U^\dagger_\nu(x-nu)
    scratch += field->component< SU3::Matrix >(Base::idx_X);          // + U_\nu(x-nu)
    field->shift(Base::idx_X, Base::dir_DOWN);

    scratch *= alpha;
    transform += scratch;
    scratch2 = localTrace(scratch);
    scratch2 *= (1.0 / 3.0);
    transform -= scratch2;

    // Dir Y
    scratch  = field->component< SU3::Matrix >(Base::idx_Y).dagger();
    scratch -= field->component< SU3::Matrix >(Base::idx_Y);

    field->shift(Base::idx_Y, Base::dir_UP);
    scratch -= field->component< SU3::Matrix >(Base::idx_Y).dagger();
    scratch += field->component< SU3::Matrix >(Base::idx_Y);
    field->shift(Base::idx_Y, Base::dir_DOWN);

    scratch *= alpha;
    transform += scratch;
    scratch2 = localTrace(scratch);
    scratch2 *= (1.0 / 3.0);
    transform -= scratch2;

    // Dir Z
    scratch  = field->component< SU3::Matrix >(Base::idx_Z).dagger();
    scratch -= field->component< SU3::Matrix >(Base::idx_Z);

    field->shift(Base::idx_Z, Base::dir_UP);
    scratch -= field->component< SU3::Matrix >(Base::idx_Z).dagger();
    scratch += field->component< SU3::Matrix >(Base::idx_Z);
    field->shift(Base::idx_Z, Base::dir_DOWN);

    scratch *= alpha;
    transform += scratch;
    scratch2 = localTrace(scratch);
    scratch2 *= (1.0 / 3.0);
    transform -= scratch2;

    Tool::reunitarize(&transform); // Possibly needs replacement!

    // Now actually transform the field
    field->component< SU3::Matrix >(Base::idx_X).leftMultiply(transform);
    transform.shift(Base::idx_X, Base::dir_DOWN);
    field->component< SU3::Matrix >(Base::idx_X).rightMultiply(transform.dagger());
    transform.shift(Base::idx_X, Base::dir_UP);

    field->component< SU3::Matrix >(Base::idx_Y).leftMultiply(transform);
    transform.shift(Base::idx_Y, Base::dir_DOWN);
    field->component< SU3::Matrix >(Base::idx_Y).rightMultiply(transform.dagger());
    transform.shift(Base::idx_Y, Base::dir_UP);

    field->component< SU3::Matrix >(Base::idx_Z).leftMultiply(transform);
    transform.shift(Base::idx_Z, Base::dir_DOWN);
    field->component< SU3::Matrix >(Base::idx_Z).rightMultiply(transform.dagger());
    transform.shift(Base::idx_Z, Base::dir_UP);

    field->component< SU3::Matrix >(Base::idx_T).leftMultiply(transform);
    transform.shift(Base::idx_T, Base::dir_DOWN);
    field->component< SU3::Matrix >(Base::idx_T).rightMultiply(transform.dagger());
    transform.shift(Base::idx_T, Base::dir_UP);

    // Calculate the current value of the gauge fixing term
    scratch  = field->component< SU3::Matrix >(Base::idx_X);
    scratch += field->component< SU3::Matrix >(Base::idx_X).dagger();
    scratch += field->component< SU3::Matrix >(Base::idx_Y);
    scratch += field->component< SU3::Matrix >(Base::idx_Y).dagger();
    scratch += field->component< SU3::Matrix >(Base::idx_Z);
    scratch += field->component< SU3::Matrix >(Base::idx_Z).dagger();

    std::complex< double > newGaugeTerm = tr(scratch);
    newGaugeTerm /= 24 * scratch.size(); // (2 * N_c * 4 * V)^-1
    if (std::abs((std::abs(newGaugeTerm) / std::abs(gaugeTerm)) - 1) < 1e-7)
      break;
    gaugeTerm = newGaugeTerm;

    // Calculate the current value of the gauge fixing term derivative (actual condition)
    // Using form from hep-lat/0104012v1
    scratch  = field->component< SU3::Matrix >(Base::idx_X);
    scratch -= field->component< SU3::Matrix >(Base::idx_X).dagger();
    scratch2  = scratch;
    scratch.shift(Base::idx_X, Base::dir_UP);
    scratch2 -= scratch;

    scratch  = field->component< SU3::Matrix >(Base::idx_Y);
    scratch -= field->component< SU3::Matrix >(Base::idx_Y).dagger();
    scratch2 += scratch;
    scratch.shift(Base::idx_Y, Base::dir_UP);
    scratch2 -= scratch;

    scratch  = field->component< SU3::Matrix >(Base::idx_Z);
    scratch -= field->component< SU3::Matrix >(Base::idx_Z).dagger();
    scratch2 += scratch;
    scratch.shift(Base::idx_Z, Base::dir_UP);
    scratch2 -= scratch;

    scratch = scratch2;
    scratch.rightMultiply(scratch2.dagger());
  }
  std::cerr << "# Finished after " << iterations << " iterations." << std::endl;
}
