
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
#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>
#include <L2/Input/FileReader.h>

#define __VINCENT__
#define __SMEARING__

int main(int argc, char **argv)
{

  size_t L = 0;
  size_t T = 0;

#ifdef __VINCENT__
 #ifdef __SMEARING__
  Input::FileReader reader("../test/my_test_input_vincent_smearing.xml");
 #else
  Input::FileReader reader("../test/my_test_input_vincent.xml");
 #endif
#else
 #ifdef __SMEARING__
  Input::FileReader reader("../test/my_test_input_smearing.xml");
 #else
  Input::FileReader reader("../test/my_test_input.xml");
 #endif
#endif

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::map< std::string, int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L, T, files, floats, positions, operators);

  Base::Weave weave(L, T);
  weave.barrier();

  std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu    = floats["mu"];
  if (weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

#ifdef __SMEARING__
  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = size_t(floats["APE_steps"]);

  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = size_t(floats["Jac_steps"]);
  if (weave.isRoot())
  {
    std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
    std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
  }
#endif

  size_t const timeslice_source = (positions[0])[Base::idx_T];
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  // size_t const source_position[4] = {0,0,0,timeslice_source};
  size_t const timeslice_boundary(T-1);
  size_t const timeslice_stochSource = (positions[1])[Base::idx_T];
  size_t const timeslice_sink = timeslice_stochSource;
  if (weave.isRoot())
    std::cout << "timeslice (stochastic wall source) = " << timeslice_stochSource << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
#ifdef __VINCENT__
 std::vector< std::string > const &smearedPropFilesDSmeared(files[2]);
 std::vector< std::string > const &smearedPropFilesUSmeared(files[3]);
#endif

//   std::vector< std::string > const &stochasticPropFilesD(files[2]);
//   std::vector< std::string > const &stochasticPropFilesU(files[3]);
//   std::vector< std::string > const &stochasticSourceFiles(files[4]);
  std::vector< std::string > const &gaugeFieldFiles(files[11]);
  std::vector< std::string > const &smearedGaugeFieldFiles(files[14]);
  std::vector< std::string > const &sourceFiles(files[12]);
  std::vector< std::string > const &smearedSourceFiles(files[13]);


  Core::Field< QCD::Gauge > gauge_field(L, T);

  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;

#ifdef __SMEARING__

  Smear::APE APE_tool(APE_alpha);
  Smear::Jacobi Jacobi_tool(Jac_alpha);

  APE_tool.smear(gauge_field, APE_iterations);


#ifdef __VINCENT__
  Tool::IO::save(&gauge_field, smearedGaugeFieldFiles[0], Tool::IO::fileILDG);
#endif

#endif

#ifdef __VINCENT__
#ifdef __SMEARING__
{
  Core::Propagator source(L, T);
  Tool::IO::load(&source, sourceFiles, Tool::IO::fileSCIDAC, 64);

  source.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  Tool::IO::save(&source, smearedSourceFiles, Tool::IO::fileSCIDAC);

  if (weave.isRoot())
  {

    std::string const inversion_command_d ("../test/invert -f ../test/invert_input_vincent_smeared_source_d >/dev/null");
    std::string const inversion_command_u ("../test/invert -f ../test/invert_input_vincent_smeared_source_u >/dev/null");
    std::string const move_command_d("mv ../../../vincent_crosscheck/files/smeared_point_sources/source.smeared.*.inverted ../../../vincent_crosscheck/files/smeared_point_source_d_propagators/");
    std::string const move_command_u("mv ../../../vincent_crosscheck/files/smeared_point_sources/source.smeared.*.inverted ../../../vincent_crosscheck/files/smeared_point_source_u_propagators/");

    int state = system(inversion_command_d.c_str());

    if(state == 0)
      std::cout << "smeared source (d) inverted successfully!" << std::endl;
    else
    {
      std::cerr << "error occurred while trying to invert smeared source (d)!" << std::endl;
      exit(1);
    }

    assert(system(move_command_d.c_str()) == 0);

    state = system(inversion_command_u.c_str());

    if(state == 0)
      std::cout << "smeared source (u) inverted successfully!" << std::endl;
    else
    {
      std::cerr << "error occurred while trying to invert smeared source (u)!" << std::endl;
      exit(1);
    }

    assert(system(move_command_u.c_str()) == 0);
  }
  weave.barrier();
}
#endif
#endif


  Core::Propagator uProp = Core::Propagator(L, T);
  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC, 64);

  if (weave.isRoot())
    std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp = Core::Propagator(L, T);
  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC, 64);

  if (weave.isRoot())
    std::cout << "d quark propagator successfully loaded\n" << std::endl;

#ifdef __VINCENT__
#ifdef __SMEARING__
  dProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  uProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
#endif
#endif



#ifdef __SMEARING__
  // Tool::IO::save(&gauge_field, gaugeFieldFiles[0] + ".smeared", Tool::IO::fileILDG);

  Core::Propagator uPropSmeared(uProp);
  Core::Propagator dPropSmeared(dProp);

  uPropSmeared.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  dPropSmeared.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  Tool::IO::save(&uPropSmeared, smearedPropFilesUSmeared, Tool::IO::fileSCIDAC);
  Tool::IO::save(&dPropSmeared, smearedPropFilesDSmeared, Tool::IO::fileSCIDAC);

#endif


//   Core::StochasticPropagator< 12 > stochastic_dProp(L, T);
//   Core::StochasticPropagator< 12 > stochastic_uProp(L, T);
//   Core::StochasticSource< 12 >     stochasticSource(L, T, Base::sou_FULLY_POLARIZED, Base::sou_PURE);
// 
//   // the load function with last parameter (size_t) precision is sort of a quick and dirty version
//   // but it also reads files without a proper header
//   Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochastic_dProp),
//                  stochasticPropFilesD, Tool::IO::fileSCIDAC, 64);
//   std::cout << "stochastic d quark propagator successfully loaded\n" << std::endl;
// 
//   Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochastic_uProp),
//                  stochasticPropFilesU, Tool::IO::fileSCIDAC, 64);
//   std::cout << "stochastic u quark propagator successfully loaded\n" << std::endl;
// 
//   Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochasticSource),
//                  stochasticSourceFiles, Tool::IO::fileSCIDAC, 64);
//   std::cout << "stochastic source successfully loaded\n" << std::endl;


  std::vector< Base::Operator > my_operators;
  my_operators.push_back(Base::op_GAMMA_4);
  // my_operators.push_back(Base::op_UNITY);

#ifdef __SMEARING__
  Core::Correlator p2p = Contract::proton_twopoint(uPropSmeared, uPropSmeared, dPropSmeared, Base::proj_NO_PROJECTOR);
#else
  Core::Correlator p2p = Contract::proton_twopoint(uProp, uProp, dProp, Base::proj_NO_PROJECTOR);
#endif
  if (weave.isRoot())
  {
    std::cout << "\nproton twopoint at t=4, full spin structure in twisted basis\n" << std::endl;
    std::cout << p2p[4] << std::endl;
    Dirac::Matrix p2pTest(p2p[4]);
    std::cout << "\nproton twopoint\n" << std::endl;
    p2p *= Base::proj_PARITY_PLUS_TM;
    std::cout << p2p << std::endl;
  }
  weave.barrier();

//   std::vector< Core::Correlator > p3p = Contract::proton_threepoint_stochastic_naive(uProp, dProp,
//                                                           stochastic_uProp, stochastic_dProp, stochasticSource,
//                                                           NULL, //no gauge field at the moment
//                                                           my_operators,
//                                                           timeslice_source, timeslice_sink);
//   p3p[0] *= Base::proj_PARITY_PLUS_TM;
//   p3p[1] *= Base::proj_PARITY_PLUS_TM;
//   std::cout.precision(8);
//   std::cout << "\nproton threepoint (stochastic naive):\n" <<std::endl;
//   std::cout << "\n d_bar*gamma0*d" <<std::endl;
//   for (size_t t=0; t<p3p[0].T(); t++)
//   {
//     if(abs(tr((p3p[0])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[0])[t])).real() << "  " << (tr((p3p[0])[t])).imag() << std::endl;
//   }
//   std::cout << "\n u_bar*gamma0*u" <<std::endl;
//   for (size_t t=0; t<p3p[1].T(); t++)
//   {
//     if(abs(tr((p3p[1])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[1])[t])).real() << "  " << (tr((p3p[1])[t])).imag() << std::endl;
//   }
//   p3p.clear();
// 
//   p3p = Contract::proton_threepoint_stochastic(uProp, dProp,
//                                                stochastic_uProp, stochastic_dProp,
//                                                stochasticSource,
//                                                timeslice_source, timeslice_stochSource,
//                                                my_operators, Base::proj_PARITY_PLUS_TM);
// 
//   std::cout.precision(8);
//   std::cout << "\nproton threepoint (stochastic):\n" <<std::endl;
//   std::cout << "\n d_bar*gamma0*d" <<std::endl;
//   for (size_t t=0; t<p3p[0].T(); t++)
//   {
//     if(abs(tr((p3p[0])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[0])[t])).real() << "  " << (tr((p3p[0])[t])).imag() << std::endl;
//   }
//   std::cout << "\n u_bar*gamma0*u" <<std::endl;
//   for (size_t t=0; t<p3p[1].T(); t++)
//   {
//     if(abs(tr((p3p[1])[t])) > 1.e-100)
//       std::cout << t << "  " << (tr((p3p[1])[t])).real() << "  " << (tr((p3p[1])[t])).imag() << std::endl;
//   }
//   p3p.clear();


#ifndef __SMEARING__
//   {
// 
//     Core::Propagator sequentialSource[16] = {
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T)};
// 
// 
//     Contract::create_sequential_source_proton_d(sequentialSource, uProp, uProp, timeslice_sink);
//     Tool::IO::save(&(sequentialSource[0]), files[9], Tool::IO::fileSCIDAC);
// 
//     Contract::create_sequential_source_proton_u(sequentialSource, dProp, uProp,timeslice_sink);
//     Tool::IO::save(&(sequentialSource[0]), files[10], Tool::IO::fileSCIDAC);
// 
//     return 0;
//   }
#endif

#ifdef __SMEARING__
//   {
// 
//     Core::Propagator sequentialSource[16] = {
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T),
//       Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T), Core::Propagator(L, T)};
// 
// 
//     Contract::create_sequential_source_proton_d(sequentialSource, uPropSmeared, uPropSmeared,
//                                                 gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
//                                                 timeslice_sink);
//     Tool::IO::save(&(sequentialSource[0]), files[9], Tool::IO::fileSCIDAC);
// 
//     Contract::create_sequential_source_proton_u(sequentialSource, dPropSmeared, uPropSmeared,
//                                                 gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
//                                                 timeslice_sink);
//     Tool::IO::save(&(sequentialSource[0]), files[10], Tool::IO::fileSCIDAC);
// 
//     //return 0;
//   }
#endif

  Core::Propagator sequentialSource(L, T);
  sequentialSource *= std::complex< double >(0, 0);

#ifdef __SMEARING__
  Contract::create_sequential_source_proton_d(sequentialSource, uPropSmeared, uPropSmeared,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_PARITY_PLUS_TM);

#else
  Contract::create_sequential_source_proton_d(sequentialSource, uProp, uProp, timeslice_sink, Base::proj_PARITY_PLUS_TM);
#endif

  Tool::IO::save(&sequentialSource, files[5], Tool::IO::fileSCIDAC);

  sequentialSource *= std::complex< double >(0, 0);

#ifdef __SMEARING__
  Contract::create_sequential_source_proton_u(sequentialSource, dPropSmeared, uPropSmeared,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_PARITY_PLUS_TM);
#else
  Contract::create_sequential_source_proton_u(sequentialSource, dProp, uProp, timeslice_sink, Base::proj_PARITY_PLUS_TM);
#endif

  Tool::IO::save(&sequentialSource, files[6], Tool::IO::fileSCIDAC);

  if (weave.isRoot())
  {

  #ifdef __VINCENT__
  #ifdef __SMEARING__
    std::string const inversion_command_d ("../test/invert -f ../test/invert_input_vincent_smeared_seq_source_d >/dev/null");
    std::string const inversion_command_u ("../test/invert -f ../test/invert_input_vincent_smeared_seq_source_u >/dev/null");
  #else
    // ...
  #endif
  #else
  #ifndef __SMEARING__
    std::string const inversion_command_d ("../test/invert -f ../test/invert_input_d >/dev/null");
    std::string const inversion_command_u ("../test/invert -f ../test/invert_input_u >/dev/null");
  #endif
  #endif

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
  }

  weave.barrier();

  Core::Propagator sequentialPropagator_u(L, T);
  Core::Propagator sequentialPropagator_d(L, T);

  Tool::IO::load(&sequentialPropagator_d, files[7], Tool::IO::fileSCIDAC);
  Tool::IO::load(&sequentialPropagator_u, files[8], Tool::IO::fileSCIDAC);

  sequentialPropagator_u.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  sequentialPropagator_d.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);

  std::vector< Core::Correlator > C3p = Contract::proton_threepoint_sequential(sequentialPropagator_u, uProp,
                                                                               sequentialPropagator_d, dProp,
                                                                               NULL, // no gauge field
                                                                               my_operators);
  if (weave.isRoot())
  {
    std::cout << "\n ubar gamma_0 u \n" << std::endl;
    std::cout << C3p[0] << std::endl;
    std::cout << "\n dbar gamma_0 d \n" << std::endl;
    std::cout << C3p[1] << std::endl;
  //   std::cout << "\n ubar u \n" << std::endl;
  //   std::cout << C3p[2] << std::endl;
  //   std::cout << "\n dbar d \n" << std::endl;
  //   std::cout << C3p[3] << std::endl;
  }

  weave.barrier();
  std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  return EXIT_SUCCESS;
}
