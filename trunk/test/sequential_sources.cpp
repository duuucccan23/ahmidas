
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

  std::cout << "testing subroutine Dirac::Matrix::outer_product(...) for coincidence of both versions" << std::endl;

  Dirac::Matrix matrixA;
  Dirac::Matrix matrixB;

  matrixA[0]=1;
  matrixA[1]=2;
  matrixA[2]=3;
  matrixA[3]=4;
  matrixA[4]=5;
  matrixA[5]=6;
  matrixA[6]=7;
  matrixA[7]=8;
  matrixA[8]=9;
  matrixA[9]=10;
  matrixA[10]=11;
  matrixA[11]=12;
  matrixA[12]=13;
  matrixA[13]=14;
  matrixA[14]=15;
  matrixA[15]=16;

  matrixB[0]=1;
  matrixB[1]=2;
  matrixB[2]=3;
  matrixB[3]=5;
  matrixB[4]=3.5;
  matrixB[5]=7;
  matrixB[6]=4.6;
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
        std::cerr << "Matrices are NOT equal in the following setup:" << std::endl;
        std::cerr << "sink   Dirac index of proton is " << alpha_f[idx_D%4] << std::endl;
        std::cerr << "source Dirac index of proton is " << alpha_i[idx_D/4] << std::endl;
        std::cerr << "Dirac::OuterProductIndexOrder is Dirac::order_FIRST_FIXED" << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  std::cout << "\nSUCCESS: Dirac::Matrix::outer_product works with Dirac::order_FIRST_FIXED\n" << std::endl;

  // result calculated by hand for alpha_f=Base:;gam_4, alpha_f=Base:;gam_4
  std::complex< double > res34[16] =  { 27.0000,  8.1000, 12.0000,  4.8000,
                                        63.0000, 18.9000, 28.0000, 11.2000,
                                        99.0000, 29.7000, 44.0000, 17.6000,
                                       135.0000, 40.5000, 60.0000, 24.0000};
  Dirac::Matrix M34(res34);

  // check whether essential subroutine works (part 2)
  {
    Dirac::OuterProductIndexOrder const order2 =  Dirac::order_OUTER_FIXED;

    matrixA.outer_product(matrixB, bufferBig, order2);

    for(size_t idx_D = 0; idx_D < 16; idx_D++)
    {
      matrixA.outer_product(matrixB, bufferSmall, order2, alpha_i[idx_D/4], alpha_f[idx_D%4]);
      if (!(Dirac::Matrix(bufferBig+(alpha_i[idx_D/4]*4+alpha_f[idx_D%4])*16) == Dirac::Matrix(bufferSmall)))
      {
        std::cerr << "\nFAILURE: matrices are NOT equal in the following setup:" << std::endl;
        std::cerr << "sink   Dirac index of proton is " << alpha_f[idx_D%4] << std::endl;
        std::cerr << "source Dirac index of proton is " << alpha_i[idx_D/4] << std::endl;
        std::cerr << "Dirac::OuterProductIndexOrder is Dirac::order_OUTER_FIXED\n" << std::endl;
        return EXIT_FAILURE;
      }
      if (alpha_f[idx_D%4] == Base::gam_3 &&  alpha_i[idx_D/4] ==  Base::gam_4)
      {
        if(!(M34 == Dirac::Matrix(bufferSmall)))
        {
          std::cerr << "\nFAILURE: control sample differs from exspected result!" << std::endl;
          std::cerr <<  "expected:\n" << M34 << "\ncalculated:\n" << Dirac::Matrix(bufferSmall) << std::endl;
          //return EXIT_FAILURE;
        }
      }
      // std::cout << "alpha_f=" << alpha_f[idx_D%4] << ", alpha_i=" << alpha_i[idx_D/4]<< std::endl;
      // std::cout << Dirac::Matrix(bufferSmall) << std::endl;
    }
  }

  std::cout << "\nSUCCESS: Dirac::Matrix::outer_product works with Dirac::order_OUTER_FIXED\n" << std::endl;


  /* ----------------------------------------------------------------*/


  double const tolerance(1.e-12);

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

  Core::Propagator uProp(L, T);
  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC);
  std::cout << "u quark propagator successfully loaded\n";

  Core::Propagator dProp(L, T);
  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC);
  std::cout << "d quark propagator successfully loaded\n";


  size_t timeslice_source(0);
  size_t timeslice_boundary(T - 1);
  uProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  dProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);

  Core::Correlator C2_P = Contract::proton_twopoint(uProp, uProp, dProp, Base::proj_NO_PROJECTOR);
  std::cout << "\nproton twopoint, full spin structure in twisted basis, at timeslice ";
  std::cout << T/2 << ":\n" << std::endl;
  Dirac::Matrix matrix_twopoint(C2_P[T/2]);
  std::cout << matrix_twopoint << std::endl;
  C2_P *= Base::proj_PARITY_PLUS_TM;
  std::cout << "\nproton twopoint, projected and traced (physical basis), at timeslice ";
  std::cout << T/2 << ":\n" << std::endl;
  std::cout << (C2_P[T/2]).trace().real() << "  " << (C2_P[T/2]).trace().imag() << std::endl;

  Core::Propagator sequentialSource[16] = {
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T)};

  Dirac::Gamma< 4 > gamma0;
  Dirac::Gamma< 5 > gamma5;

  // this is a good test for the sequential source generation
  // just multiplying the sequential source (without gamma_5 and dagger) to a propagator gives the twopoint

  Core::Propagator sequentialSource_fixedProjector_d(L, T);
  Core::Propagator sequentialSource_fixedProjector_u(L, T);


  Core::Propagator dProp_mod(dProp);
  Core::Propagator::iterator it = dProp_mod.begin();
  while(it != dProp_mod.end())
  {
    (*it).right_multiply_proton();
    (*it).left_multiply_proton();
    (*it).transposeFull();
    ++it;
  }

  sequentialSource_fixedProjector_d *= std::complex< double >(0, 0);
  sequentialSource_fixedProjector_u *= std::complex< double >(0, 0);
  Base::Weave weave(L, T);

  QCD::Tensor tmp[16];

  size_t localIndex;
  for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
  {
    for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
    {
      for(size_t idx_X = 0; idx_X < L; idx_X++)
      {
        localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, T/2);
        /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
        if (localIndex == weave.localVolume())
          continue;
        QCD::make_sequential_d(tmp, uProp[localIndex], uProp[localIndex]);
        QCD::make_sequential_d(sequentialSource_fixedProjector_d[localIndex], uProp[localIndex], uProp[localIndex], Base::proj_PARITY_PLUS_TM);
        QCD::make_sequential_u(sequentialSource_fixedProjector_u[localIndex], dProp_mod[localIndex], uProp[localIndex], Base::proj_PARITY_PLUS_TM);
        for(size_t idx_D = 0; idx_D < 16; idx_D++)
          (sequentialSource[idx_D])[localIndex] = tmp[idx_D];
      }
    }
  }

  std::cout.precision(6);

  Dirac::Matrix matrix_u, matrix_d;

  for(size_t idx_D = 0; idx_D < 16; idx_D++)
  {
    Core::Correlator p2p_seqD(L, T, (sequentialSource[idx_D]).contract(dProp_mod));
    //Core::Correlator p2p_seqU(L, T, (sequentialSource[idx_D]).contract(uProp));
    p2p_seqD.sumOverSpatialVolume();
    //p2p_seqU.sumOverSpatialVolume();
    matrix_d[idx_D] = (p2p_seqD[2]).trace();
    //matrix_u[idx_D] = (p2p_seqU[2]).trace();
  }
  std::cout << "\nproton twopoint from sequential source (d), full spin structure in twisted basis\n" << std::endl;
  std::cout << matrix_d << std::endl;
  Dirac::Matrix matrixTest_d (matrix_d);
  std::cout << "\nthis is the trace of the projected and traced twopoint in physical basis\n" << std::endl;
  Dirac::Matrix matrix2(gamma5*matrix_d);
  matrix2 *= std::complex< double >(0, 1);
  matrix_d = gamma0*matrix_d;
  matrix_d += matrix2;
  matrix_d *= 0.5;
  std::cout.width(16);
  std::cout << matrix_d.trace().real();
  std::cout.width(16);
  std::cout << matrix_d.trace().imag();
  std::cout << std::endl;
  std::cout << "\nproton twopoint from sequential source (d), fixed projector\n" << std::endl;
  Core::Correlator p2p_seq(L, T, sequentialSource_fixedProjector_d.contract(dProp_mod));
  p2p_seq.sumOverSpatialVolume();
  std::cout << p2p_seq << "\n" << std::endl;

  if(abs(matrix_d.trace() - p2p_seq[T/2].trace())/abs(matrix_d.trace()) < tolerance)
  {
    std::cout << "\nSUCCESS: sequential method with free proton Dirac indices gives the same result as method with fixed ones.\n" << std::endl;
  }
  else
  {
    std::cerr << "\nFAILURE: sequential method with free proton Dirac indices gives result DIFFERENT from method with fixed ones.\n" << std::endl;
    return EXIT_FAILURE;
  }

  if(abs(p2p_seq[T/2].trace() - C2_P[T/2].trace())/abs(C2_P[T/2].trace()) < tolerance)
  {
    std::cout << "\nSUCCESS: twopoint at particular timeslice can be recovered from sequential source (d) times d quark propagator\n" << std::endl;
  }
  else
  {
    std::cerr << "\nFAILURE: twopoint at particular timeslice can NOT be recovered from sequential source (d) times d quark propagator\n" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "\nproton twopoint from sequential source (u), fixed projector\n" << std::endl;
  Core::Correlator p2p_seq_u(L, T, sequentialSource_fixedProjector_u.contract(dProp));
  p2p_seq_u.sumOverSpatialVolume();
  std::cout << p2p_seq_u << "\n" << std::endl;


  if(abs((p2p_seq_u[T/2]).trace() - (C2_P[T/2]).trace())/abs(C2_P[T/2].trace()) < tolerance)
  {
    std::cout << "\nSUCCESS: twopoint at particular timeslice can be recovered from sequential source (u) times u quark propagator\n" << std::endl;
  }
  else
  {
    std::cerr << "\nFAILURE: twopoint at particular timeslice can NOT be recovered from sequential source (u) times u quark propagator\n" << std::endl;
    return EXIT_FAILURE;
  }


  return EXIT_FAILURE;
}