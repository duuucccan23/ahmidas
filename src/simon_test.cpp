
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
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>
#include <L2/Input/FileReader.h>


int main(int argc, char **argv)
{

  size_t L = 0;
  size_t T = 0;

  Input::FileReader reader("../test/my_test_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::map< std::string, int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L, T, files, floats, positions, operators);

  std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu = floats["mu"];

  std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

  size_t const timeslice_source = (positions[0])[Base::idx_T];
  std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  // size_t const source_position[4] = {0,0,0,timeslice_source};
  size_t const timeslice_boundary(T-1);
  size_t const timeslice_stochSource = (positions[1])[Base::idx_T];
  size_t const timeslice_sink = timeslice_stochSource;
  std::cout << "timeslice (stochastic wall source) = " << timeslice_stochSource << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &stochasticPropFilesD(files[2]);
  std::vector< std::string > const &stochasticPropFilesU(files[3]);
  std::vector< std::string > const &stochasticSourceFiles(files[4]);


//   std::cout << "The following files are going to be read:" << std::endl;
// 
//   for (int f=0; f<12; f++)
//   {
//       std::cout << propfilesU[f]           << "\n" << propfilesD[f]           << std::endl;
//       std::cout << stochasticPropFilesD[f] << "\n" << stochasticPropFilesU[f] << std::endl;
//   }
// 
//   for (int f=0; f<12; f++)
//   {
//     std::cout << stochasticSourceFiles[f] << std::endl;
//   }

  Core::Propagator uProp = Core::Propagator(L, T);

  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC, 64);

  std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp = Core::Propagator(L, T);

  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC, 64);

  std::cout << "d quark propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 12 > stochastic_dProp(L, T);
  Core::StochasticPropagator< 12 > stochastic_uProp(L, T);
  Core::StochasticSource< 12 >     stochasticSource(L, T, Base::sou_FULLY_POLARIZED, Base::sou_PURE);

  // the load function with last parameter (size_t) precision is sort of a quick and dirty version
  // but it also reads files without a proper header
  Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochastic_dProp),
                 stochasticPropFilesD, Tool::IO::fileSCIDAC, 64);
  std::cout << "stochastic d quark propagator successfully loaded\n" << std::endl;

  Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochastic_uProp),
                 stochasticPropFilesU, Tool::IO::fileSCIDAC, 64);
  std::cout << "stochastic u quark propagator successfully loaded\n" << std::endl;

  Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochasticSource),
                 stochasticSourceFiles, Tool::IO::fileSCIDAC, 64);
  std::cout << "stochastic source successfully loaded\n" << std::endl;


  std::vector< Base::Operator > my_operators;
  my_operators.push_back(Base::op_GAMMA_4);
  // my_operators.push_back(Base::op_UNITY);

  Core::Correlator p2p = Contract::proton_twopoint(uProp, uProp, dProp, Base::proj_NO_PROJECTOR);
  std::cout << "\nproton twopoint at t=4, full spin structure in twisted basis\n" << std::endl;
  std::cout << p2p[4] << std::endl;
  Dirac::Matrix p2pTest(p2p[4]);
  std::cout << "\nproton twopoint\n" << std::endl;
  p2p *= Base::proj_PARITY_PLUS_TM;
  std::cout << p2p << std::endl;


  std::vector< Core::Correlator > p3p = Contract::proton_threepoint_stochastic_naive(uProp, dProp,
                                                          stochastic_uProp, stochastic_dProp, stochasticSource,
                                                          NULL, //no gauge field at the moment
                                                          my_operators,
                                                          timeslice_source, timeslice_sink);
  p3p[0] *= Base::proj_PARITY_PLUS_TM;
  p3p[1] *= Base::proj_PARITY_PLUS_TM;
  std::cout.precision(8);
  std::cout << "\nproton threepoint (stochastic naive):\n" <<std::endl;
  std::cout << "\n d_bar*gamma0*d" <<std::endl;
  for (size_t t=0; t<p3p[0].T(); t++)
  {
    if(abs(tr((p3p[0])[t])) > 1.e-100)
      std::cout << t << "  " << (tr((p3p[0])[t])).real() << "  " << (tr((p3p[0])[t])).imag() << std::endl;
  }
  std::cout << "\n u_bar*gamma0*u" <<std::endl;
  for (size_t t=0; t<p3p[1].T(); t++)
  {
    if(abs(tr((p3p[1])[t])) > 1.e-100)
      std::cout << t << "  " << (tr((p3p[1])[t])).real() << "  " << (tr((p3p[1])[t])).imag() << std::endl;
  }
  p3p.clear();

  p3p = Contract::proton_threepoint_stochastic(uProp, dProp,
                                               stochastic_uProp, stochastic_dProp,
                                               stochasticSource,
                                               timeslice_source, timeslice_stochSource,
                                               my_operators, Base::proj_PARITY_PLUS_TM);

  std::cout.precision(8);
  std::cout << "\nproton threepoint (stochastic):\n" <<std::endl;
  std::cout << "\n d_bar*gamma0*d" <<std::endl;
  for (size_t t=0; t<p3p[0].T(); t++)
  {
    if(abs(tr((p3p[0])[t])) > 1.e-100)
      std::cout << t << "  " << (tr((p3p[0])[t])).real() << "  " << (tr((p3p[0])[t])).imag() << std::endl;
  }
  std::cout << "\n u_bar*gamma0*u" <<std::endl;
  for (size_t t=0; t<p3p[1].T(); t++)
  {
    if(abs(tr((p3p[1])[t])) > 1.e-100)
      std::cout << t << "  " << (tr((p3p[1])[t])).real() << "  " << (tr((p3p[1])[t])).imag() << std::endl;
  }
  p3p.clear();

//   std::cout << "\n d_bar*1*d" <<std::endl;
//   for (size_t t=0; t<p3p[2].T(); t++)
//   {
//     if(abs(tr((p3p[2])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[2])[t])).real() << "  " << (tr((p3p[2])[t])).imag() << std::endl;
//   }
//   std::cout << "\n u_bar*1*u" <<std::endl;
//   for (size_t t=0; t<p3p[3].T(); t++)
//   {
//     if(abs(tr((p3p[3])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[3])[t])).real() << "  " << (tr((p3p[3])[t])).imag() << std::endl;
//   }
/*
  p3p.clear();*/

  {
    Core::Propagator sequentialSource[16] = {
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T)};

    // Contract::create_sequential_source_proton_d(sequentialSource, uProp, uProp, timeslice_sink);
    // Tool::IO::save(&(sequentialSource[0]), files[9], Tool::IO::fileSCIDAC);

    //std::cerr << "sequential source (d)\n\n" << sequentialSource[0] << std::endl;

    Dirac::Gamma< 4 > gamma0;
    Dirac::Gamma< 5 > gamma5;

  // this is a good test for the sequential source generation
  // just multiplying the sequential source (without gamma_5 and dagger) to a propagator gives the twopoint

    Core::Propagator sequentialSource_fixedProjector(L, T);

    Core::Propagator dProp_mod(dProp);
    {
      Core::Propagator::iterator it = dProp_mod.begin();
      while(it != dProp_mod.end())
      {
        (*it).right_multiply_proton();
        (*it).left_multiply_proton();
        (*it).transposeFull();
        ++it;
      }
    }

    sequentialSource_fixedProjector *= std::complex< double >(0, 0);
    Base::Weave weave(L, T);

    QCD::Tensor tmp[16];

    size_t localIndex;
    for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L; idx_X++)
        {
          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, timeslice_sink);
          /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
          if (localIndex == weave.localVolume())
            continue;
          QCD::make_sequential_d(tmp, uProp[localIndex], uProp[localIndex]);
          QCD::make_sequential_d(sequentialSource_fixedProjector[localIndex], uProp[localIndex], uProp[localIndex], Base::proj_PARITY_PLUS_TM);
          for(size_t idx_D = 0; idx_D < 16; idx_D++)
            (sequentialSource[idx_D])[localIndex] = tmp[idx_D];
        }
      }
    }

    std::cout.precision(6);

    Dirac::Matrix matrix1;

    for(size_t idx_D = 0; idx_D < 16; idx_D++)
    {
      Core::Correlator p2p_seq(L, T, (sequentialSource[idx_D]).contract(dProp_mod));
      p2p_seq.sumOverSpatialVolume();
      matrix1[idx_D] = (p2p_seq[4]).trace();
    }
    std::cout << "\nproton twopoint from sequential source (d), full spin structure in twisted basis\n" << std::endl;
    std::cout << matrix1 << std::endl;
    Dirac::Matrix matrixTest (matrix1);
    std::cout << "\nthis is the trace of the projected and traced twopoint in physical basis\n" << std::endl;
    Dirac::Matrix matrix2(gamma5*matrix1);
    matrix2 *= std::complex< double >(0, 1);
    matrix1 = gamma0*matrix1;
    matrix1 += matrix2;
    matrix1 *= 0.5;
    std::cout.width(16);
    std::cout << matrix1.trace().real();
    std::cout.width(16);
    std::cout << matrix1.trace().imag();
    std::cout << std::endl;
    std::cout << "\nproton twopoint from sequential source (d), fixed projector\n" << std::endl;
    Core::Correlator p2p_seq(L, T, sequentialSource_fixedProjector.contract(dProp_mod));
    p2p_seq.sumOverSpatialVolume();
    std::cout << p2p_seq << "\n" << std::endl;

    matrixTest -= p2pTest;
    for(size_t idx = 0; idx < 16; idx++)
      matrixTest[idx] = abs(matrixTest[idx])/abs(p2pTest[idx]);

    std::cout << "\nrelative difference between correct and this solution\n" << std::endl;
    std::cout << matrixTest << std::endl;

    // now the same for the u current sequential source

    sequentialSource_fixedProjector *= std::complex< double >(0, 0);

    for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L; idx_X++)
        {
          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, timeslice_sink);
          /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
          if (localIndex == weave.localVolume())
            continue;
          QCD::make_sequential_u(sequentialSource_fixedProjector[localIndex], dProp_mod[localIndex], uProp[localIndex], Base::proj_PARITY_PLUS_TM);
          QCD::make_sequential_u(tmp, dProp_mod[localIndex], uProp[localIndex]);
          for(size_t idx_D = 0; idx_D < 16; idx_D++)
          {
            (sequentialSource[idx_D])[localIndex] = tmp[idx_D];
          }
        }
      }
    }

    for(size_t idx_D = 0; idx_D < 16; idx_D++)
    {
      Core::Correlator p2p_seq(L, T, (sequentialSource[idx_D]).contract(uProp));
      p2p_seq.sumOverSpatialVolume();
      matrix1[idx_D] = (p2p_seq[4]).trace()*0.5;
    }
    std::cout << "\nproton twopoint from sequential source (u), full spin structure in twisted basis\n" << std::endl;
    std::cout << matrix1 << std::endl;
    matrixTest = matrix1;
    std::cout << "\nthis is the trace of the projected and traced twopoint in physical basis\n" << std::endl;
    matrix2= gamma5*matrix1;
    matrix2 *= std::complex< double >(0, 1);
    matrix1 = gamma0*matrix1;
    matrix1 += matrix2;
    matrix1 *= 0.5;
    std::cout.width(16);
    std::cout << matrix1.trace().real();
    std::cout.width(16);
    std::cout << matrix1.trace().imag();
    std::cout << std::endl;
    std::cout << "\nproton twopoint from sequential source (u), fixed projector\n" << std::endl;
    Core::Correlator p2p_seq_u(L, T, sequentialSource_fixedProjector.contract(uProp));
    p2p_seq_u.sumOverSpatialVolume();
    p2p_seq_u *= 0.5;
    std::cout << p2p_seq_u << "\n" << std::endl;

    matrixTest -= p2pTest;
    for(size_t idx = 0; idx < 16; idx++)
      matrixTest[idx] = abs(matrixTest[idx])/abs(p2pTest[idx]);

    std::cout << "\nrelative difference between correct and this solution\n" << std::endl;
    std::cout << matrixTest << std::endl;

  }

  {

    Core::Propagator sequentialSource[16] = {
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
      Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T)};

    Core::Propagator dProp_mod(dProp);
    Core::Propagator::iterator it = dProp_mod.begin();
    while(it != dProp_mod.end())
    {
      (*it).right_multiply_proton();
      (*it).left_multiply_proton();
      (*it).transposeFull();
      ++it;
    }

    Contract::create_sequential_source_proton_d(sequentialSource, uProp, uProp, timeslice_sink);
    Tool::IO::save(&(sequentialSource[0]), files[9], Tool::IO::fileSCIDAC);

    Contract::create_sequential_source_proton_u(sequentialSource, dProp, uProp, timeslice_sink);
    Tool::IO::save(&(sequentialSource[0]), files[10], Tool::IO::fileSCIDAC);
  }

  Core::Propagator sequentialSource(L, T);
  sequentialSource *= std::complex< double >(0, 0);

  Contract::create_sequential_source_proton_d(sequentialSource, uProp, uProp, timeslice_sink, Base::proj_PARITY_PLUS_TM_STAR);

  Tool::IO::save(&sequentialSource, files[5], Tool::IO::fileSCIDAC);

  Contract::create_sequential_source_proton_u(sequentialSource, uProp, dProp, timeslice_sink, Base::proj_PARITY_PLUS_TM_STAR);

  Tool::IO::save(&sequentialSource, files[6], Tool::IO::fileSCIDAC);

  std::string const inversion_command_d ("../test/invert -f ../test/invert_input_d >/dev/null");
  std::string const inversion_command_u ("../test/invert -f ../test/invert_input_u >/dev/null");

  int state = system(inversion_command_d.c_str());

  if(state == 0)
    std::cout << "sequential source (d) inverted successfully!" << std::endl;
  else
  {
    std::cerr << "error occurred while trying to invert sequential source (d)!" << std::endl;
    exit(1);
  }

  state = system(inversion_command_u.c_str());

  if(state == 0)
    std::cout << "sequential source (u) inverted successfully!" << std::endl;
  else
  {
    std::cerr << "error occurred while trying to invert sequential source (u)!" << std::endl;
    exit(1);
  }

  Core::Propagator sequentialPropagator_u(L, T);
  Core::Propagator sequentialPropagator_d(L, T);

  Tool::IO::load(&sequentialPropagator_d, files[7], Tool::IO::fileSCIDAC);
  Tool::IO::load(&sequentialPropagator_u, files[8], Tool::IO::fileSCIDAC);

  sequentialPropagator_u.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  sequentialPropagator_d.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);

  std::vector< Core::Correlator > C3p = Contract::proton_threepoint_sequential(sequentialPropagator_u, uProp,
                                                                 sequentialPropagator_d, dProp,
                                                                 my_operators);
  std::cout << "\n ubar gamma_0 u \n" << std::endl;
  std::cout << C3p[0] << std::endl;
  std::cout << "\n dbar gamma_0 d \n" << std::endl;
  std::cout << C3p[1] << std::endl;
//   std::cout << "\n ubar u \n" << std::endl;
//   std::cout << C3p[2] << std::endl;
//   std::cout << "\n dbar d \n" << std::endl;
//   std::cout << C3p[3] << std::endl;

  std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  return EXIT_SUCCESS;
}
