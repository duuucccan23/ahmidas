#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

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
#define __GAMMA_1__
// use operator involving gamma_2
// #define __GAMMA_2__
// use operator involving gamma_3
// #define __GAMMA_3__
// use operator involving gamma_4
// #define __GAMMA_4__

int main(int argc, char **argv)
{

  size_t L = 0;
  size_t T = 0;


  Input::FileReader reader("./generate_sequential_sources_alternative_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::map< std::string, int > operators;
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


  size_t const sourceSinkSeparation = size_t(floats["t_source_insertion"]);
  assert(sourceSinkSeparation > 0 && sourceSinkSeparation < T-1);


  size_t const * const source_position = positions[0];
  size_t const timeslice_source = source_position[Base::idx_T] % T;
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  size_t const t_op = (timeslice_source +  sourceSinkSeparation) % T;
  if (weave.isRoot())
    std::cout << "timeslice (operator insertion) = " << t_op << std::endl;

  // make sure the boundary is not crossed by source-sink correlaton function
  size_t const timeslice_boundary = (timeslice_source + (T/2)) % T;;
  if (weave.isRoot())
    std::cout << "timeslice (boundary) = " << timeslice_boundary << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
  std::vector< std::string > const &seqSourceFilesD(files[3]);
  std::vector< std::string > const &seqSourceFilesU(files[4]);

  Core::Field< QCD::Gauge > gauge_field(L, T);

  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;


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


  {

  Core::Propagator sequentialSource_u(L, T);
  Core::Propagator sequentialSource_d(L, T);


#ifdef __GAMMA_1__
  Contract::create_sequential_source_proton_fixed_insertion_timeslice(&sequentialSource_u, &sequentialSource_d,
                                                                        *uProp, *dProp, t_op, Base::op_GAMMA_15);
#endif  

#ifdef __GAMMA_2__
  Contract::create_sequential_source_proton_fixed_insertion_timeslice(&sequentialSource_u, &sequentialSource_d,
                                                                        *uProp, *dProp, t_op, Base::op_GAMMA_25);
#endif  

#ifdef __GAMMA_3__
  Contract::create_sequential_source_proton_fixed_insertion_timeslice(&sequentialSource_u, &sequentialSource_d,
                                                                        *uProp, *dProp, t_op, Base::op_GAMMA_35);
#endif  

#ifdef __GAMMA_4__
  Contract::create_sequential_source_proton_fixed_insertion_timeslice(&sequentialSource_u, &sequentialSource_d,
                                                                        *uProp, *dProp, t_op, Base::op_GAMMA_45);
#endif  

  Tool::IO::save(&sequentialSource_d, seqSourceFilesD, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource_u, seqSourceFilesU, Tool::IO::fileSCIDAC);

  delete dProp;
  delete uProp;

  double norm_u(sequentialSource_u.norm());

  if (weave.isRoot())
  {
    std::cout.precision(10);
    std::cout << std::scientific << "norm of sequential source (u): " << norm_u << std::endl;
  }

  double norm_d(sequentialSource_d.norm());

  if (weave.isRoot())
  {
    std::cout.precision(10);
    std::cout << std::scientific << "norm of sequential source (d): " << norm_d << std::endl;
  }

  if (weave.isRoot())
    std::cout << "sequential sources generated and saved successfully\n" << std::endl;

  }
  return EXIT_SUCCESS;
}
