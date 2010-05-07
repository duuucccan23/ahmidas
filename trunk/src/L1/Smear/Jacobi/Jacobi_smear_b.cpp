#include "Jacobi.ih"

namespace Smear
{
  void Jacobi::smear(Core::Field< QCD::Tensor > *tensorField, Core::Field< QCD::Gauge > &gaugeField) const
  {
    Core::Field< QCD::Tensor > neighbour(*tensorField);
    neighbour *= d_kappa;

    // Sum the field with the parallelly transported neighbours.
    (*tensorField) += Transport::step(neighbour, gaugeField, Base::idx_X, Base::dir_DOWN);
    (*tensorField) += Transport::step(neighbour, gaugeField, Base::idx_Y, Base::dir_DOWN);
    (*tensorField) += Transport::step(neighbour, gaugeField, Base::idx_Z, Base::dir_DOWN);

    (*tensorField) += Transport::step(neighbour, gaugeField, Base::idx_X, Base::dir_UP);
    (*tensorField) += Transport::step(neighbour, gaugeField, Base::idx_Y, Base::dir_UP);
    (*tensorField) += Transport::step(neighbour, gaugeField, Base::idx_Z, Base::dir_UP);

    (*tensorField) *= d_weight;
  }
}
