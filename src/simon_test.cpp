
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

#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__

int main(int argc, char **argv)
{

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("../test/my_test_input_fixed_timeslice.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::map< std::string, int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);

  size_t const L(L_tmp);
  size_t const T(T_tmp);

  Base::Weave weave(L, T);
  weave.barrier();

  std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu    = floats["mu"];
  if (weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = size_t(floats["APE_steps"]);

  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = size_t(floats["Jac_steps"]);
  if (weave.isRoot())
  {
    std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
    std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
  }

  size_t const timeslice_source = (positions[0])[Base::idx_T];
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  // size_t const source_position[4] = {0,0,0,timeslice_source};
  size_t const timeslice_boundary(T-1);
  size_t const t_op = size_t(floats["t_insertion"]);
  if (weave.isRoot())
    std::cout << "timeslice (operator insertion) = " << t_op << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
  std::vector< std::string > const &sequentialSourceFilesD(files[3]);
  std::vector< std::string > const &sequentialSourceFilesU(files[4]);
  std::vector< std::string > const &sequentialPropFilesD(files[5]);
  std::vector< std::string > const &sequentialPropFilesU(files[6]);

  Core::Field< QCD::Gauge > gauge_field(L, T);

  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;


  if (APE_iterations > 0)
  {
    Smear::APE APE_tool(APE_alpha);
    APE_tool.smear(gauge_field, APE_iterations);
  }

  Smear::Jacobi Jacobi_tool(Jac_alpha);


  Core::Propagator uProp(L, T);
  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC, 64);

  if (weave.isRoot())
    std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp(L, T);
  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC, 64);

  if (weave.isRoot())
    std::cout << "d quark propagator successfully loaded\n" << std::endl;

#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
  dProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  uProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
#endif

  Core::Correlator p2p = Contract::proton_twopoint(uProp, uProp, dProp, Base::proj_NO_PROJECTOR);

  if (weave.isRoot())
  {
    std::cout << "\nproton twopoint\n" << std::endl;
    p2p *= Base::proj_PARITY_PLUS_TM;
    std::cout << p2p << std::endl;
  }
  weave.barrier();

  {

  Core::Propagator sequentialSource_u(L, T);
  Core::Propagator sequentialSource_d(L, T);


  Contract::create_sequential_source_proton_fixed_insertion_timeslice(&sequentialSource_u, &sequentialSource_d,
                                                            uProp, dProp, t_op, Base::op_GAMMA_4);


  Tool::IO::save(&sequentialSource_d, sequentialSourceFilesD, Tool::IO::fileSCIDAC);
  Tool::IO::save(&sequentialSource_u, sequentialSourceFilesU, Tool::IO::fileSCIDAC);

  }

  if (weave.isRoot())
  {

    //std::string const inversion_command_d ("../test/invert -f ../test/simon/invert_input_d >/dev/null");
    //std::string const inversion_command_u ("../test/invert -f ../test/simon/invert_input_u >/dev/null");
    std::string const inversion_command_d ("../test/invert -f ../test/simon/invert_input_d");
    std::string const inversion_command_u ("../test/invert -f ../test/simon/invert_input_u");

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

  Tool::IO::load(&sequentialPropagator_d, sequentialPropFilesD, Tool::IO::fileSCIDAC);
  Tool::IO::load(&sequentialPropagator_u, sequentialPropFilesU, Tool::IO::fileSCIDAC);

#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
  sequentialPropagator_u.changeBoundaryConditions_uniformToFixed(t_op, timeslice_boundary);
  sequentialPropagator_d.changeBoundaryConditions_uniformToFixed(t_op, timeslice_boundary);
#endif


  if (Jac_iterations > 0)
  {
    uProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
    dProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
    sequentialPropagator_u.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
    sequentialPropagator_d.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  }

  Core::Correlator C3p_d  = Contract::proton_twopoint(uProp, uProp, sequentialPropagator_d, Base::proj_PARITY_PLUS_TM);
  Core::Correlator C3p_u1 = Contract::proton_twopoint(uProp, sequentialPropagator_u, dProp, Base::proj_PARITY_PLUS_TM);
  Core::Correlator C3p_u2 = Contract::proton_twopoint(sequentialPropagator_u, uProp, dProp, Base::proj_PARITY_PLUS_TM);

  C3p_u1 += C3p_u2;

  if (weave.isRoot())
  {
    std::cout << "\n ubar gamma_0 u \n" << std::endl;
    std::cout << C3p_u1 << std::endl;
    std::cout << "\n dbar gamma_0 d \n" << std::endl;
    std::cout << C3p_d << std::endl;
  }

  weave.barrier();
  std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  return EXIT_SUCCESS;
}
