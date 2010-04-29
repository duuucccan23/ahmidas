
#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Base/Weave.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Matrix.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>


int main(int argc, char **argv)
{


  Dirac::Matrix matrixA;
  Dirac::Matrix matrixB;

  matrixA[0]=1;
  matrixA[1]=2;
  matrixA[2]=std::complex<double> (3, 1);
  matrixA[3]=4;
  matrixA[4]=5;
  matrixA[5]=6;
  matrixA[6]=7;
  matrixA[7]=8;
  matrixA[8]=std::complex<double> (9, 2);
  matrixA[9]=10;
  matrixA[10]=11;
  matrixA[11]=12;
  matrixA[12]=13;
  matrixA[13]=std::complex<double> (14, 3);
  matrixA[14]=15;
  matrixA[15]=16;

  matrixB[0]=std::complex<double> (1, 4);
  matrixB[1]=2;
  matrixB[2]=3;
  matrixB[3]=5;
  matrixB[4]=3.5;
  matrixB[5]=7;
  matrixB[6]=std::complex<double> (2, 3);
  matrixB[7]=5.1;
  matrixB[8]=3.9;
  matrixB[9]=11;
  matrixB[10]=1.2;
  matrixB[11]=2.3;
  matrixB[12]=9;
  matrixB[13]=2.7;
  matrixB[14]=4;
  matrixB[15]=1.6;

  std::vector < Base::DiracIndex > alpha_i;
  alpha_i.push_back(Base::gam_1);
  alpha_i.push_back(Base::gam_2);
  alpha_i.push_back(Base::gam_3);
  alpha_i.push_back(Base::gam_4);
  std::vector < Base::DiracIndex > alpha_f;
  alpha_f.push_back(Base::gam_1);
  alpha_f.push_back(Base::gam_2);
  alpha_f.push_back(Base::gam_3);
  alpha_f.push_back(Base::gam_4);

  std::complex< double > bufferSmall[16];
  std::complex< double > bufferBig[16*16];

  // check whether essential subroutine works (part 1)
  {
    Dirac::OuterProductIndexOrder const order1 =  Dirac::order_FIRST_FIXED;

    matrixA.outer_product(matrixB, bufferBig, order1);

    for(size_t idx_D = 0; idx_D < 16; idx_D++)
    {
      matrixA.outer_product(matrixB, bufferSmall, order1, alpha_i[idx_D/4], alpha_f[idx_D%4]);
      if (!(Dirac::Matrix(bufferBig+(alpha_i[idx_D/4]*4+alpha_f[idx_D%4])*16) == Dirac::Matrix(bufferSmall)))
      {
        std::cout << "Matrices are NOT equal" << std::endl;
        std::cout << "sink Dirac index of proton is" << alpha_f[idx_D%4] << std::endl;
        std::cout << "sink Dirac index of proton is" << alpha_i[idx_D/4] << std::endl;
        std::cout << "Dirac::OuterProductIndexOrder is" << order1 << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  std::cout << "Dirac::Matrix::outer_product works with Dirac::order_FIRST_FIXED" << std::endl;

  // check whether essential subroutine works (part 2)
  {
    Dirac::OuterProductIndexOrder const order2 =  Dirac::order_OUTER_FIXED;

    matrixA.outer_product(matrixB, bufferBig, order2);

    for(size_t idx_D = 0; idx_D < 16; idx_D++)
    {
      matrixA.outer_product(matrixB, bufferSmall, order2, alpha_i[idx_D/4], alpha_f[idx_D%4]);
      if (!(Dirac::Matrix(bufferBig+(alpha_i[idx_D/4]*4+alpha_f[idx_D%4])*16) == Dirac::Matrix(bufferSmall)))
      {
        std::cout << "Matrices are NOT equal" << std::endl;
        std::cout << "sink Dirac index of proton is" << alpha_f[idx_D%4] << std::endl;
        std::cout << "sink Dirac index of proton is" << alpha_i[idx_D/4] << std::endl;
        std::cout << "Dirac::OuterProductIndexOrder is" << order2 << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  std::cout << "Dirac::Matrix::outer_product works with Dirac::order_OUTER_FIXED" << std::endl;


  /* ----------------------------------------------------------------*/


  const size_t L = 4;
  const size_t T = 4;

  std::vector<std::string> propfilesU;
  std::vector<std::string> propfilesD;

  const std::string filename_base1("../../test/source4x4");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss <<  ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    propfilesU.push_back(std::string(filename_base1).append("_u").append(oss.str()));
    propfilesD.push_back(std::string(filename_base1).append("_d").append(oss.str()));
  }

  Core::Propagator *uProp = new Core::Propagator(L, T);
  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);
  std::cout << "u quark propagator successfully loaded\n";

  Core::Propagator *dProp = new Core::Propagator(L, T);
  Tool::IO::load(dProp, propfilesD, Tool::IO::fileSCIDAC);
  std::cout << "d quark propagator successfully loaded\n";


  size_t timeslice_source(0);
  size_t timeslice_boundary(T - 1);
  uProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  dProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);

  Core::Correlator C2_P = Contract::proton_twopoint(*uProp, *uProp, *dProp, Base::proj_PARITY_PLUS_TM);




  return EXIT_FAILURE;
}