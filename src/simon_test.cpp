
#include <cstring>
#include <vector>
#include <map>
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
  // size_t const timeslice_boundary(T-1);
  size_t const timeslice_stochSource = (positions[1])[Base::idx_T];
  std::cout << "timeslice (stochastic wall source) = " << timeslice_stochSource << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
//   std::vector< std::string > const &stochasticPropFilesD(files[2]);
//   std::vector< std::string > const &stochasticPropFilesU(files[3]);
//   std::vector< std::string > const &stochasticSourceFiles(files[4]);


  std::cout << "The following files are going to be read:" << std::endl;

  for (int f=0; f<12; f++)
  {
      std::cout << propfilesU[f]           << "\n" << propfilesD[f]           << std::endl;
//       std::cout << stochasticPropFilesD[f] << "\n" << stochasticPropFilesU[f] << std::endl;
  }

  for (int f=0; f<12; f++)
  {
//     std::cout << stochasticSourceFiles[f] << std::endl;
  }

  Core::Propagator uProp = Core::Propagator(L, T);

  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC, 64);

  std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp = Core::Propagator(L, T);

  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC, 64);


  std::cout << "d quark propagator successfully loaded\n" << std::endl;


//   Core::StochasticPropagator< 12 > *stochastic_dProp = new Core::StochasticPropagator< 12 >(L, T);
//   Core::StochasticPropagator< 12 > *stochastic_uProp = new Core::StochasticPropagator< 12 >(L, T);
//   Core::StochasticSource< 12 > *stochasticSource = new Core::StochasticSource< 12 >(L, T, Base::sou_FULLY_POLARIZED,
//                                                                                           Base::sou_PURE);
//   // the load function with last parameter (size_t) precision is sort of a quick and dirty version
//   // but it also reads files without a proper header
//   Tool::IO::load(dynamic_cast< Core::Propagator *> (stochastic_dProp),
//                  stochasticPropFilesD, Tool::IO::fileSCIDAC, 64);
//   std::cout << "stochastic d quark propagator successfully loaded\n" << std::endl;
// 
// 
//   Tool::IO::load(dynamic_cast< Core::Propagator *> (stochastic_uProp),
//                  stochasticPropFilesU, Tool::IO::fileSCIDAC, 64);
//   std::cout << "stochastic u quark propagator successfully loaded\n" << std::endl;
// 
// 
//   Tool::IO::load(dynamic_cast< Core::Propagator *> (stochasticSource),
//                  stochasticSourceFiles, Tool::IO::fileSCIDAC, 64);
//   std::cout << "stochastic source successfully loaded\n" << std::endl;


//   std::vector< Base::Operator > my_operators;
//   my_operators.push_back(Base::op_GAMMA_4);
// 
//   std::vector< Core::Correlator > p3p = Contract::proton_threepoint_stochastic(*uProp, *dProp,
//                                                                      *stochastic_uProp, *stochastic_dProp,
//                                                                      *stochasticSource,
//                                                                      timeslice_source, timeslice_stochSource,
//                                                                      my_operators, Base::proj_PARITY_PLUS_TM);
//   std::cout << "\nproton threepoint:\n" <<std::endl;
//   std::cout << "\n d_bar*Op*d" <<std::endl;
//   for (size_t t=0; t<p3p[0].T(); t++)
//   {
//    // if(abs(tr((p3p[0])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[0])[t])).real() << "  " << (tr((p3p[0])[t])).imag() << std::endl;
//   }
//   std::cout << "\n u_bar*Op*u" <<std::endl;
//   for (size_t t=0; t<p3p[1].T(); t++)
//   {
//     //if(abs(tr((p3p[1])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[1])[t])).real() << "  " << (tr((p3p[1])[t])).imag() << std::endl;
//   }


//   std::cout << "\nthat is supposed to be the result:\n"
//     << "0   1.16145514e-03  -7.64689531e-05\n"
//     << "1   1.37746719e-03  -1.02786839e-05\n"
//     << "2   4.29020380e-05   3.89420319e-08\n"
//     << "3   1.34970449e-03   1.41088319e-05"
//     << std::endl;


  Core::Propagator sequentialSource[16] = {Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T),
                                            Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T),
                                            Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T),
                                            Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T), Core:: Propagator(L,T)};

  Contract::create_sequential_source_proton_d(sequentialSource, uProp, uProp);


  std::cout << sequentialSource[0] << std::endl;


  //Core::Correlator C3pd = Contract::proton_threepoint_d(Core:: Propagator * const bw_prop, Core:: Propagator * const fw_prop);

  std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

// #ifdef __MPI_ARCH__
//   if (!MPI::Is_finalized())
//     MPI::Finalize();
// #endif

  return EXIT_SUCCESS;
}
