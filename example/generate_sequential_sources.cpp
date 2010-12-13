#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Ahmidas.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>
#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>
#include <L2/Input/FileReader.h>


// comment this in if propagators have uniform temporal boundary contitions
// (e.g. the HMC inverter does this)
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS__

// use operator involving gamma_1
// #define __GAMMA_1__
// use operator involving gamma_2
// #define __GAMMA_2__
// use operator involving gamma_3
// #define __GAMMA_3__
// use operator involving gamma_4
 #define __GAMMA_4__

// use all projectors
//#define __ALL_GAMMA__

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  size_t L = 0;
  size_t T = 0;


  Input::FileReader reader("./generate_sequential_sources_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::vector< int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L, T, files, floats, positions, operators);

  Base::Weave weave(L, T);
  weave.barrier();

  if (weave.isRoot())
    std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu    = floats["mu"];
  if (weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;


  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = size_t(floats["APE_steps"]);
  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = size_t(floats["Jac_steps"]);

  size_t const sourceSinkSeparation = size_t(floats["sourceSinkSeparation"]);
  assert(sourceSinkSeparation > 0 && sourceSinkSeparation < T-1);

  if (weave.isRoot())
  {
    std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
    std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
  }

  size_t const * const source_position = positions[0];
  size_t const timeslice_source = source_position[Base::idx_T]  % T;
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  size_t const timeslice_sink = (timeslice_source +  sourceSinkSeparation) % T;
  if (weave.isRoot())
    std::cout << "timeslice (sink) = " << timeslice_sink << std::endl;

  // make sure the boundary is not crossed by source-sink correlaton function
  size_t const timeslice_boundary = (timeslice_source + (T/2)) % T;;
  if (weave.isRoot())
    std::cout << "timeslice (boundary) = " << timeslice_boundary << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
#ifdef __ALL_GAMMA__
  std::vector< std::string > const &seqSourceFilesD1(files[3]);
  std::vector< std::string > const &seqSourceFilesU1(files[4]);
  std::vector< std::string > const &seqSourceFilesD2(files[5]);
  std::vector< std::string > const &seqSourceFilesU2(files[6]);
  std::vector< std::string > const &seqSourceFilesD3(files[7]);
  std::vector< std::string > const &seqSourceFilesU3(files[8]);
  std::vector< std::string > const &seqSourceFilesD4(files[9]);
  std::vector< std::string > const &seqSourceFilesU4(files[10]);
#else  
  std::vector< std::string > const &seqSourceFilesD(files[3]);
  std::vector< std::string > const &seqSourceFilesU(files[4]);
#endif


  Core::Field< QCD::Gauge > gauge_field(L, T);

  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;


  Smear::APE APE_tool(APE_alpha);
  Smear::Jacobi Jacobi_tool(Jac_alpha);

  APE_tool.smear(gauge_field, APE_iterations);


  Core::Propagator *uProp = new Core::Propagator(L, T);
  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

  if (weave.isRoot())
    std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator *dProp = new Core::Propagator(L, T);
  Tool::IO::load(dProp, propfilesD, Tool::IO::fileSCIDAC);


  if (weave.isRoot())
    std::cout << "d quark propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS__
  dProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  uProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
#endif


  uProp->smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  dProp->smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  if (weave.isRoot())
    std::cout << "propagators smeared successfully\n" << std::endl;


#ifdef __ALL_GAMMA__
  Core::Propagator sequentialSource1(L, T);
  sequentialSource1 *= std::complex< double >(0, 0); // initialize with zero
  Core::Propagator sequentialSource2(L, T);
  sequentialSource2 *= std::complex< double >(0, 0); // initialize with zero
  Core::Propagator sequentialSource3(L, T);
  sequentialSource3 *= std::complex< double >(0, 0); // initialize with zero
  Core::Propagator sequentialSource4(L, T);
  sequentialSource4 *= std::complex< double >(0, 0); // initialize with zero
  
  Contract::create_sequential_source_proton_u(sequentialSource1, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_1_PLUS_TM);
  Contract::create_sequential_source_proton_u(sequentialSource2, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_2_PLUS_TM);
  Contract::create_sequential_source_proton_u(sequentialSource3, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_3_PLUS_TM);
  Contract::create_sequential_source_proton_u(sequentialSource4, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_PARITY_PLUS_TM);
#else
  
  Core::Propagator sequentialSource(L, T);
  sequentialSource *= std::complex< double >(0, 0); // initialize with zero
  
#ifdef __GAMMA_1__
  Contract::create_sequential_source_proton_u(sequentialSource, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_1_PLUS_TM);
#endif

#ifdef __GAMMA_2__
  Contract::create_sequential_source_proton_u(sequentialSource, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_2_PLUS_TM);
#endif

#ifdef __GAMMA_3__
  Contract::create_sequential_source_proton_u(sequentialSource, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_3_PLUS_TM);
#endif

#ifdef __GAMMA_4__
  Contract::create_sequential_source_proton_u(sequentialSource, *dProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_PARITY_PLUS_TM);
#endif
#endif
  
  delete dProp;

#ifdef __ALL_GAMMA__
  Tool::IO::save(&sequentialSource1, seqSourceFilesU1, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource2, seqSourceFilesU2, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource3, seqSourceFilesU3, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource4, seqSourceFilesU4, Tool::IO::fileSCIDAC);
#else
  Tool::IO::save(&sequentialSource, seqSourceFilesU, Tool::IO::fileSCIDAC);
#endif

#ifdef __ALL_GAMMA__
  double norm_u1(sequentialSource1.norm());
  double norm_u2(sequentialSource2.norm());
  double norm_u3(sequentialSource3.norm());
  double norm_u4(sequentialSource4.norm());
#else
  double norm_u(sequentialSource.norm());
#endif

  if (weave.isRoot())
  {
    std::cout.precision(10);
#ifdef __ALL_GAMMA__
    std::cout << std::scientific << "norm of sequential source (u, projector 1): " << norm_u1 << std::endl;
    std::cout << std::scientific << "norm of sequential source (u, projector 2): " << norm_u2 << std::endl;
    std::cout << std::scientific << "norm of sequential source (u, projector 3): " << norm_u3 << std::endl;
    std::cout << std::scientific << "norm of sequential source (u, projector 4): " << norm_u4 << std::endl;
#else
    std::cout << std::scientific << "norm of sequential source (u): " << norm_u << std::endl;
#endif    
  }

#ifdef __ALL_GAMMA__

  sequentialSource1 *= std::complex< double >(0, 0); // initialize with zero
  sequentialSource2 *= std::complex< double >(0, 0); // initialize with zero
  sequentialSource3 *= std::complex< double >(0, 0); // initialize with zero
  sequentialSource4 *= std::complex< double >(0, 0); // initialize with zero
  
  Contract::create_sequential_source_proton_d(sequentialSource1, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_1_PLUS_TM);
  Contract::create_sequential_source_proton_d(sequentialSource2, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_2_PLUS_TM);
  Contract::create_sequential_source_proton_d(sequentialSource3, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_3_PLUS_TM);
  Contract::create_sequential_source_proton_d(sequentialSource4, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_PARITY_PLUS_TM);
#else

  sequentialSource *= std::complex< double >(0, 0); // initialize with zero

#ifdef __GAMMA_1__
  Contract::create_sequential_source_proton_d(sequentialSource, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_1_PLUS_TM);
#endif

#ifdef __GAMMA_2__
  Contract::create_sequential_source_proton_d(sequentialSource, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_2_PLUS_TM);
#endif

#ifdef __GAMMA_3__
  Contract::create_sequential_source_proton_d(sequentialSource, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_3_PLUS_TM);
#endif

#ifdef __GAMMA_4__
  Contract::create_sequential_source_proton_d(sequentialSource, *uProp, *uProp,
                                              gauge_field, Smear::sm_Jacobi, Jac_iterations, Jac_alpha,
                                              timeslice_sink, Base::proj_PARITY_PLUS_TM);
#endif

#endif

  delete uProp;

#ifdef __ALL_GAMMA__
  Tool::IO::save(&sequentialSource1, seqSourceFilesD1, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource2, seqSourceFilesD2, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource3, seqSourceFilesD3, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource4, seqSourceFilesD4, Tool::IO::fileSCIDAC);
#else
  Tool::IO::save(&sequentialSource, seqSourceFilesD, Tool::IO::fileSCIDAC);
#endif

#ifdef __ALL_GAMMA__
  double norm_d1(sequentialSource1.norm());
  double norm_d2(sequentialSource2.norm());
  double norm_d3(sequentialSource3.norm());
  double norm_d4(sequentialSource4.norm());
#else
  double norm_d(sequentialSource.norm());
#endif

  if (weave.isRoot())
  {
    std::cout.precision(10);
#ifdef __ALL_GAMMA__
    std::cout << std::scientific << "norm of sequential source (d, projector 1): " << norm_d1 << std::endl;
    std::cout << std::scientific << "norm of sequential source (d, projector 2): " << norm_d2 << std::endl;
    std::cout << std::scientific << "norm of sequential source (d, projector 3): " << norm_d3 << std::endl;
    std::cout << std::scientific << "norm of sequential source (d, projector 4): " << norm_d4 << std::endl;
#else
    std::cout << std::scientific << "norm of sequential source (d): " << norm_d << std::endl;
#endif    
  }

  if (weave.isRoot())
    std::cout << "sequential sources generated and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
