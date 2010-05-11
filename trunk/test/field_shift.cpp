
#include <complex>
#include <iomanip>
#include <iostream>

#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>
#include <L0/Core/Field.h>

int main(int argc, char **argv)
{

  size_t errors(0);

  size_t const L =  8;
  size_t const T = 16;

  Base::Weave weave(L, T);

  Core::Field< QCD::Gauge > gauge_field(L, T);

  size_t idx(0);

  SU3::Matrix links[4];

  size_t localIndex;
  for(size_t idx_T = 0; idx_T < T; idx_T++)
  {
    for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L; idx_X++)
        {
          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);

          // globalCoordToLocalIndex returns local volume if local data is not available on this cpu
          if (localIndex == weave.localVolume())
          {
            ++idx;
            continue;
          }

          links[0] = SU3::Matrix(std::complex< double >(idx, 0));
          links[1] = SU3::Matrix(std::complex< double >(idx, 1));
          links[2] = SU3::Matrix(std::complex< double >(idx, 2));
          links[3] = SU3::Matrix(std::complex< double >(idx, 3));

          gauge_field[localIndex] = QCD::Gauge(links);

          ++idx;
        }
      }
    }
  }


  idx = 0;
  gauge_field.shift(Base::idx_T, Base::dir_DOWN);
  gauge_field.shift(Base::idx_T, Base::dir_DOWN);
  gauge_field.shift(Base::idx_Y, Base::dir_UP);
  gauge_field.shift(Base::idx_X, Base::dir_UP);
  gauge_field.shift(Base::idx_Z, Base::dir_DOWN);
  gauge_field.shift(Base::idx_T, Base::dir_UP);
  gauge_field.shift(Base::idx_Y, Base::dir_DOWN);
  gauge_field.shift(Base::idx_T, Base::dir_UP);
  gauge_field.shift(Base::idx_X, Base::dir_DOWN);
  gauge_field.shift(Base::idx_Z, Base::dir_UP);

  for(size_t idx_T = 0; idx_T < T; idx_T++)
  {
    for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L; idx_X++)
        {

          weave.barrier();

          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);

          // globalCoordToLocalIndex returns local volume if local data is not available on this cpu
          if (localIndex == weave.localVolume())
          {
            ++idx;
            continue;
          }

          std::cout.width(4);
          std::cout << idx << " -> ";
          std::cout.width(4);
          std::cout  << (((gauge_field[localIndex])[0])(0, 0)).real() << std::endl;

          if (!(((gauge_field[localIndex])[0])(0, 0) == std::complex< double >(idx, 0)
            && ((gauge_field[localIndex])[1])(0, 0) == std::complex< double >(idx, 1)
            && ((gauge_field[localIndex])[2])(0, 0) == std::complex< double >(idx, 2)
            && ((gauge_field[localIndex])[3])(0, 0) == std::complex< double >(idx, 3)))
          {
            ++errors;
          }

          ++idx;
        }
      }
    }
  }



  if (errors == 0)
  {
    std::cout << "SUCCESS: shifting along some closed path preserves the field!" << std::endl;
    return EXIT_SUCCESS;
  }

  std::cerr << "FAILURE: shifting along some closed path results a different field!" << std::endl;
  return EXIT_FAILURE;
}
