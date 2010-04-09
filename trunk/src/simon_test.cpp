
#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Base/Knuth.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>

// #define __MPI_ARCH__
// #define __STOCHASTIC_TWOPOINT__

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
  std::cout << "timeslice (stochastic wall source) = " << timeslice_stochSource << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &stochasticPropFilesD(files[2]);
  std::vector< std::string > const &stochasticPropFilesU(files[3]);
  std::vector< std::string > const &stochasticSourceFiles(files[4]);


  std::cout << "The following files are going to be read:" << std::endl;

  for (int f=0; f<12; f++)
  {
      std::cout << propfilesU[f]           << "\n" << propfilesD[f]           << std::endl;
      std::cout << stochasticPropFilesD[f] << "\n" << stochasticPropFilesU[f] << std::endl;
  }

  for (int f=0; f<12; f++)
  {
    std::cout << stochasticSourceFiles[f] << std::endl;
  }


  /* prepare random gauge transformation */
//   int seed = 83141764;
//   Base::Knuth::instance(seed);
//   Core::Field< SU3::Matrix > randomGaugeTrafo(L, T);
//   Core::Field< SU3::Matrix >::iterator I_rgt = randomGaugeTrafo.begin();
//   while (I_rgt != randomGaugeTrafo.end())
//   {
//     (*I_rgt).setToRandom();
//     //std::cout << *I_rgt << std::endl;
//     ++I_rgt;
//   }

  /* function to be called:
     gaugeTransform_fixedSource(randomGaugeTrafo, source_position)
  */


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


//   // perform random gauge transformation of point source propagators
//   dProp->gaugeTransform_fixedSource(randomGaugeTrafo, source_position);
//   uProp->gaugeTransform_fixedSource(randomGaugeTrafo, source_position);


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

//     std::cout << *stochasticSource << std::endl;

//   uProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
//   dProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);


  // perform field shift ...
  // stochastic_dProp->shift(Base::idx_T, Base::dir_UP, int(timeslice_stochSource)-int(timeslice_source));

#ifdef __STOCHASTIC_TWOPOINT__

  Core::Propagator *dProp_stoch =  new Core::Propagator((stochasticSource->createStochasticPropagator_fixedSink(*stochastic_uProp, source_position)).revert());

  Core::Propagator *uProp_stoch =  new Core::Propagator((stochasticSource->createStochasticPropagator_fixedSink(*stochastic_dProp, source_position)).revert());

//   std::ofstream f("fake_propagator");
//   f << *dProp_stoch << std::endl;
//   f.close();


//   delete stochasticSource;
//   delete stochastic_dProp;


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
  delete dProp_stoch;
  delete uProp_stoch;
#endif

  std::vector< Base::Operator > my_operators;
  my_operators.push_back(Base::op_GAMMA_4);

  std::vector< Core::Correlator > p3p = Contract::proton_threepoint_stochastic(*uProp, *dProp,
                                                                     *stochastic_uProp, *stochastic_dProp,
                                                                     *stochasticSource,
                                                                     timeslice_source, timeslice_stochSource,
                                                                     my_operators, Base::proj_PARITY_PLUS_TM);
  std::cout << "\nproton threepoint:\n" <<std::endl;
  std::cout << "\n d_bar*Op*d" <<std::endl;
  for (size_t t=0; t<p3p[0].getT(); t++)
  {
   // if(abs(tr((p3p[0])[t])) > 1.e-100)
      std::cout << t << "  " << (tr((p3p[0])[t])).real() << "  " << (tr((p3p[0])[t])).imag() << std::endl;
  }
  std::cout << "\n u_bar*Op*u" <<std::endl;
  for (size_t t=0; t<p3p[1].getT(); t++)
  {
    //if(abs(tr((p3p[1])[t])) > 1.e-100)
      std::cout << t << "  " << (tr((p3p[1])[t])).real() << "  " << (tr((p3p[1])[t])).imag() << std::endl;
  }


//   std::cout << "\nthat is supposed to be the result:\n"
//     << "0   1.16145514e-03  -7.64689531e-05\n"
//     << "1   1.37746719e-03  -1.02786839e-05\n"
//     << "2   4.29020380e-05   3.89420319e-08\n"
//     << "3   1.34970449e-03   1.41088319e-05"
//     << std::endl;

  delete uProp;
  delete dProp;


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
