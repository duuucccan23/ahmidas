#include <L1/Smear/APE/APE.ih>

namespace Smear
{
  void APE::smear(Core::Field< QCD::Gauge > &field) const
  {
    Core::Field< SU3::Matrix > xAPE = Path::staple(field, Base::idx_Y, Base::dir_UP, Base::idx_X, Base::dir_UP);
    xAPE += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xAPE += Path::staple(field, Base::idx_Z, Base::dir_UP, Base::idx_X, Base::dir_UP);
    xAPE += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xAPE *= d_alpha;

    Core::Field< SU3::Matrix > yAPE = Path::staple(field, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP);
    yAPE += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yAPE += Path::staple(field, Base::idx_Z, Base::dir_UP, Base::idx_Y, Base::dir_UP);
    yAPE += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yAPE *= d_alpha;

    Core::Field< SU3::Matrix > zAPE = Path::staple(field, Base::idx_X, Base::dir_UP, Base::idx_Z, Base::dir_UP);
    zAPE += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zAPE += Path::staple(field, Base::idx_Y, Base::dir_UP, Base::idx_Z, Base::dir_UP);
    zAPE += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zAPE *= d_alpha;

    field.component< SU3::Matrix >(Base::idx_X) += xAPE;
    field.component< SU3::Matrix >(Base::idx_Y) += yAPE;
    field.component< SU3::Matrix >(Base::idx_Z) += zAPE;

    Tool::reunitarize(&field);
  }
}
