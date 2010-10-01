#include "Tensor.ih"

namespace QCD
{

  void Tensor::setToRandom(Base::SourcePolarization const DState, Base::SourceColorState const CState,
                           Base::SourceStochasticTypeFlag const type, uint64_t seed)
  {
    std::fill_n(d_data, 144, std::complex< double >(0, 0));
    std::complex< double > tmp_data [12];
    if (seed == 0)
    {
      switch (type)
      {
        case Base::sou_Z4:
        {
          std::generate_n(reinterpret_cast< double* >(tmp_data), 24, Base::Random::Z2);
          double const norm = 1.0/sqrt(2.0);
          std::transform(reinterpret_cast< double* >(tmp_data), reinterpret_cast< double* >(tmp_data) + 24,
                        reinterpret_cast< double* >(tmp_data), std::bind2nd(std::multiplies< double >(), norm));
          break;
        }
        case Base::sou_Z2:
        {
          std::generate_n(tmp_data, 12, Base::Random::Z2);
  //         std::generate_n(reinterpret_cast< double* >(tmp_data), 24, Base::Random::Z2);
  //         for (size_t idx=1; idx<24, idx+=2)
  //           (reinterpret_cast< double* >(tmp_data))[idx] = 0.0; // set imaginary parts to zero
          break;
        }
        case Base::sou_P1:
        {
          std::fill_n(tmp_data, 12, std::complex< double >(1, 0));
          break;
        }
        case Base::sou_M1:
        {
          std::fill_n(tmp_data, 12, std::complex< double >(-1, 0));
          break;
        }
        default:
          std::cerr << "This Base::SourceColorState is not implemented in Tensor::setToRandom_Z4()\n";
          std::cerr << "Aborting ..." << std::endl;
          exit(1);
      }
    }
    else
    {
      switch (type)
      {
        case Base::sou_Z4:
        {
          std::generate_n(reinterpret_cast< double* >(tmp_data), 24, Base::Z2(1.0, seed));
          double const norm = 1.0/sqrt(2.0);
          std::transform(reinterpret_cast< double* >(tmp_data), reinterpret_cast< double* >(tmp_data) + 24,
                        reinterpret_cast< double* >(tmp_data), std::bind2nd(std::multiplies< double >(), norm));
          break;
        }
        case Base::sou_Z2:
        {
          std::generate_n(tmp_data, 12, Base::Z2(1.0, seed));
  //         std::generate_n(reinterpret_cast< double* >(tmp_data), 24, Base::Random::Z2);
  //         for (size_t idx=1; idx<24, idx+=2)
  //           (reinterpret_cast< double* >(tmp_data))[idx] = 0.0; // set imaginary parts to zero
          break;
        }
        case Base::sou_P1:
        {
          std::fill_n(tmp_data, 12, std::complex< double >(1, 0));
          break;
        }
        case Base::sou_M1:
        {
          std::fill_n(tmp_data, 12, std::complex< double >(-1, 0));
          break;
        }
        default:
          std::cerr << "This Base::SourceColorState is not implemented in Tensor::setToRandom_Z4()\n";
          std::cerr << "Aborting ..." << std::endl;
          exit(1);
      }
    }
    switch (DState)
    {
      // case Base::sou_UNPOLARIZED:
      case Base::sou_PARTLY_POLARIZED:
        switch (CState)
        {
          // in this case we have no dilution at all
          case Base::sou_GENERIC:
            std::copy(tmp_data, tmp_data + 12, d_data);
            break;
          default:
          std::cerr << "This Base::SourceColorState is not implemented in Tensor::setToRandom_Z4()\n";
          std::cerr << "Aborting ..." << std::endl;
          exit(1);
        }
        break;
      case Base::sou_FULLY_POLARIZED:
        switch (CState)
        {
          // this is what is called spin diluted
          case Base::sou_GENERIC:
            std::copy(tmp_data,     tmp_data +  3, d_data);
            std::copy(tmp_data + 3, tmp_data +  6, d_data + 39);
            std::copy(tmp_data + 6, tmp_data +  9, d_data + 78);
            std::copy(tmp_data + 9, tmp_data + 12, d_data + 117);
            break;
          // this is what is called spin and color diluted
          case Base::sou_PURE:
            d_data[  0] = tmp_data[ 0];
            d_data[ 13] = tmp_data[ 1];
            d_data[ 26] = tmp_data[ 2];
            d_data[ 39] = tmp_data[ 3];
            d_data[ 52] = tmp_data[ 4];
            d_data[ 65] = tmp_data[ 5];
            d_data[ 78] = tmp_data[ 6];
            d_data[ 91] = tmp_data[ 7];
            d_data[104] = tmp_data[ 8];
            d_data[117] = tmp_data[ 9];
            d_data[130] = tmp_data[10];
            d_data[143] = tmp_data[11];
            break;
          default:
          std::cerr << "This Base::SourceColorState is not implemented in Tensor::setToRandom_Z4()\n";
          std::cerr << "Aborting ..." << std::endl;
          exit(1);
        }
        break;
      default:
      std::cerr << "This Base::SourcePolarization is not implemented in Tensor::setToRandom_Z4()\n";
      std::cerr << "Aborting ..." << std::endl;
      exit(1);
    }
  }
}
