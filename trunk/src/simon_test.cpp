
#include <cstring>
#include <vector>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>

// #define __MPI_ARCH__


int main(int argc, char **argv)
{

// #ifdef __MPI_ARCH__
//   MPI::Init(argc, argv);
//   int numprocs(MPI::COMM_WORLD.Get_size());
//   int myid(MPI::COMM_WORLD.Get_rank());
// #endif

  const size_t L = 4;
  const size_t T = 8;

  std::vector<std::string> propfilesU;
  std::vector<std::string> propfilesD;
  std::vector<std::string> stochasticPropFilesD;
  std::vector<std::string> stochasticPropFilesU;
  std::vector<std::string> stochasticSourceFiles;

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
  std::cout << "The following files are going to be read:" << std::endl;

  const std::string filename_baseU("../../../crosscheck_dru/point_source_u_propagators/source");
  const std::string filename_baseD("../../../crosscheck_dru/point_source_d_propagators/source");
  // this is necessary since we have the stochastic source at the sink, reverting propagator changes flavour
  const std::string filename_baseA("../../../crosscheck_dru/stochastic_source_u_propagators/source.0000.0000");
  const std::string filename_baseB("../../../crosscheck_dru/stochastic_source_d_propagators/source.0000.0000");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss <<  ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    propfilesU.push_back(std::string(filename_baseU).append(oss.str()));
    propfilesD.push_back(std::string(filename_baseD).append(oss.str()));
    stochasticPropFilesD.push_back(std::string(filename_baseA).append(oss.str()));
    stochasticPropFilesU.push_back(std::string(filename_baseB).append(oss.str()));
#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << propfilesU[f]           << "\n" << propfilesD[f]           << std::endl;
      std::cout << stochasticPropFilesD[f] << "\n" << stochasticPropFilesU[f] << std::endl;
  }

  const std::string filename_base2("../../../crosscheck_dru/stochastic_sources/source.0000");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss <<  ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss.flush();
    stochasticSourceFiles.push_back(std::string(filename_base2).append(oss.str()));
#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << stochasticSourceFiles[f] << std::endl;
  }


  Core::Propagator *uProp = new Core::Propagator(L, T);

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC, 64);

  #ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator *dProp = new Core::Propagator(L, T);

  Tool::IO::load(dProp, propfilesD, Tool::IO::fileSCIDAC, 64);

  #ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << "d quark propagator successfully loaded\n" << std::endl;


  Core::StochasticPropagator< 12 > *stochastic_dProp = new Core::StochasticPropagator< 12 >(L, T);
  Core::StochasticPropagator< 12 > *stochastic_uProp = new Core::StochasticPropagator< 12 >(L, T);
  Core::StochasticSource< 12 > *stochasticSource = new Core::StochasticSource< 12 >(L, T, Base::sou_FULLY_POLARIZED,
                                                                                          Base::sou_PURE);
  // the load function with last parameter (size_t) precision is sort of a quick and dirty version
  // but it also reads files without a proper header
  Tool::IO::load(dynamic_cast< Core::Propagator *> (stochastic_dProp),
                 stochasticPropFilesD, Tool::IO::fileSCIDAC, 64);
#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
    std::cout << "stochastic d quark propagator successfully loaded\n" << std::endl;


  Tool::IO::load(dynamic_cast< Core::Propagator *> (stochastic_uProp),
                 stochasticPropFilesU, Tool::IO::fileSCIDAC, 64);
#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
    std::cout << "stochastic u quark propagator successfully loaded\n" << std::endl;


  Tool::IO::load(dynamic_cast< Core::Propagator *> (stochasticSource),
                 stochasticSourceFiles, Tool::IO::fileSCIDAC, 64);
#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
    std::cout << "stochastic source successfully loaded\n" << std::endl;



  size_t const timeslice_source(0);
  size_t const source_position[4] = {0,0,0,timeslice_source};
//   size_t timeslice_boundary(T-1);
  size_t const timeslice_stochSource(4);
//   uProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
//   dProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);


  // perform field shift ...
  // stochastic_dProp->shift(Base::idx_T, Base::dir_UP, int(timeslice_stochSource)-int(timeslice_source));

  Core::Propagator *dProp_stoch =  new Core::Propagator((stochasticSource->createStochasticPropagator_fixedSink(*stochastic_dProp, source_position)).revert());

  Core::Propagator *uProp_stoch =  new Core::Propagator((stochasticSource->createStochasticPropagator_fixedSink(*stochastic_uProp, source_position)).revert());

//   std::ofstream f("fake_propagator");
//   f << *dProp_stoch << std::endl;
//   f.close();


  delete stochasticSource;
  delete stochastic_dProp;


  Core::Correlator C2_P = Contract::proton_twopoint(*uProp, *uProp, *dProp, Base::proj_PARITY_PLUS_TM);

  Core::Correlator C2_P_stoch_d  = Contract::proton_twopoint(*uProp,       *uProp,       *dProp_stoch, Base::proj_PARITY_PLUS_TM);
  Core::Correlator C2_P_stoch_u1 = Contract::proton_twopoint(*uProp_stoch, *uProp,       *dProp,       Base::proj_PARITY_PLUS_TM);
  Core::Correlator C2_P_stoch_u2 = Contract::proton_twopoint(*uProp,       *uProp_stoch, *dProp,       Base::proj_PARITY_PLUS_TM);

  std::cout << "\nstandard proton twopoint:\n" <<std::endl;
  std::ofstream fout("p2p.dat");
  for (size_t t=0; t<C2_P.getT(); t++)
  {
    fout      << std::scientific << std::setprecision(8) << std::showpos;
    std::cout << std::scientific << std::setprecision(8) << std::showpos;
    fout      << t << " " << (tr(C2_P[t])).real() << " " << (tr(C2_P[t])).imag() << std::endl;
    std::cout << t << "  " << (tr(C2_P[t])).real() << "  " << (tr(C2_P[t])).imag() << std::endl;
    //std::cout << C2_P[t] << std::endl;
  }

  std::cout << "\nproton twopoint using stochastic d line:\n" <<std::endl;
  for (size_t t=0; t<C2_P_stoch_d.getT(); t++)
  {
    if(abs(tr(C2_P_stoch_d[t])) > 1.e-100)
     std::cout << t << "  " << (tr(C2_P_stoch_d[t])).real() << "  " << (tr(C2_P_stoch_d[t])).imag() << std::endl;
  }
  std::cout << "\nproton twopoint using stochastic u line (1):\n" <<std::endl;
  for (size_t t=0; t<C2_P_stoch_u1.getT(); t++)
  {
    if(abs(tr(C2_P_stoch_u1[t])) > 1.e-100)
      std::cout << t << "  " << (tr(C2_P_stoch_u1[t])).real() << "  " << (tr(C2_P_stoch_u1[t])).imag() << std::endl;
  }
  std::cout << "\nproton twopoint using stochastic u line (2):\n" <<std::endl;
  for (size_t t=0; t<C2_P_stoch_u2.getT(); t++)
  {
    if(abs(tr(C2_P_stoch_u2[t])) > 1.e-100)
      std::cout << t << "  " << (tr(C2_P_stoch_u2[t])).real() << "  " << (tr(C2_P_stoch_u2[t])).imag() << std::endl;
  }



  fout.close();

//   std::cout << "\nthat is supposed to be the result:\n"
//     << "0   1.16145514e-03  -7.64689531e-05\n"
//     << "1   1.37746719e-03  -1.02786839e-05\n"
//     << "2   4.29020380e-05   3.89420319e-08\n"
//     << "3   1.34970449e-03   1.41088319e-05"
//     << std::endl;

  delete uProp;
  delete dProp;
  delete dProp_stoch;

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
  std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

// #ifdef __MPI_ARCH__
//   if (!MPI::Is_finalized())
//     MPI::Finalize();
// #endif

  return EXIT_SUCCESS;
}
