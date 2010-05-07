#include "Jacobi.ih"

namespace Smear
{
  void Jacobi::smear(Core::Field< QCD::Spinor > *spinorField, Core::Field< QCD::Gauge > &gaugeField) const
  {
    Core::Field< QCD::Spinor > neighbourSpinor(*spinorField);
    neighbourSpinor *= d_kappa;

    // Sum the field with the parallelly transported neighbours.
    (*spinorField) += Transport::step(neighbourSpinor, gaugeField, Base::idx_X, Base::dir_DOWN);
    (*spinorField) += Transport::step(neighbourSpinor, gaugeField, Base::idx_Y, Base::dir_DOWN);
    (*spinorField) += Transport::step(neighbourSpinor, gaugeField, Base::idx_Z, Base::dir_DOWN);

    (*spinorField) += Transport::step(neighbourSpinor, gaugeField, Base::idx_X, Base::dir_UP);
    (*spinorField) += Transport::step(neighbourSpinor, gaugeField, Base::idx_Y, Base::dir_UP);
    (*spinorField) += Transport::step(neighbourSpinor, gaugeField, Base::idx_Z, Base::dir_UP);

    (*spinorField) *= d_weight;
  }
}
