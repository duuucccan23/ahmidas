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


// comment this if you don't want the twopoint function to be calculated
#define __CALCULATE_TWOPOINT__

int main(int argc, char **argv)
{

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("../example/contract_3point_sequential_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::map< std::string, int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);

  size_t const L(L_tmp);
  size_t const T(T_tmp);

  Base::Weave weave(L, T);

  if (weave.isRoot());
    std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu    = floats["mu"];
  if (weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

  size_t const sourceSinkSeparation = size_t(floats["sourceSinkSeparation"]);
  assert(sourceSinkSeparation > 0 && sourceSinkSeparation < T-1);


  size_t const * const source_position = positions[0];
  size_t const timeslice_source = source_position[Base::idx_T];
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  size_t const timeslice_sink = (timeslice_source +  sourceSinkSeparation) % T;
  if (weave.isRoot())
    std::cout << "timeslice (sink) = " << timeslice_sink << std::endl;

  // make sure the boundary is not crossed by source-sink correlaton function
  size_t const timeslice_boundary = timeslice_sink > timeslice_source ? (T-1) : timeslice_source-1;
  if (weave.isRoot())
    std::cout << "timeslice (boundary) = " << timeslice_boundary << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
  std::vector< std::string > const &seqPropFilesD(files[3]);
  std::vector< std::string > const &seqPropFilesU(files[4]);

  Core::Field< QCD::Gauge > gauge_field(L, T);
  if (weave.isRoot());
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;

  Core::Propagator forwardProp_u(L, T);
  Tool::IO::load(&forwardProp_u, propfilesU, Tool::IO::fileSCIDAC, 64);
  if (weave.isRoot());
    std::cout << "u quark forward propagator successfully loaded\n" << std::endl;

  Core::Propagator forwardProp_d(L, T);
  Tool::IO::load(&forwardProp_d, propfilesD, Tool::IO::fileSCIDAC, 64);
  if (weave.isRoot())
    std::cout << "d quark forward propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS__
//   forwardProp_d.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
//   forwardProp_u.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
#endif

  Core::Propagator backwardProp_u(L, T);
  Core::Propagator backwardProp_d(L, T);

  Tool::IO::load(&backwardProp_u, seqPropFilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark backward propagator successfully loaded\n" << std::endl;

  Tool::IO::load(&backwardProp_d, seqPropFilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark backward propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS__
  backwardProp_u.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
#endif


  std::vector< Base::Operator > my_operators;

  my_operators.push_back(Base::op_GAMMA_4);
  my_operators.push_back(Base::op_O44);
  my_operators.push_back(Base::op_O11);
  my_operators.push_back(Base::op_O22);
  my_operators.push_back(Base::op_O33);


  std::vector< Core::Correlator > C3p = Contract::proton_threepoint_sequential(backwardProp_u, forwardProp_u,
                                                                               backwardProp_d, forwardProp_d,
                                                                               &gauge_field,
                                                                               my_operators);

  if (weave.isRoot())
  {
    std::ofstream fout("output_3point.dat");
    std::cout << "\n ubar gamma_0 u \n" << std::endl;
    std::cout << C3p[0] << std::endl;
    fout << C3p[0] << std::endl;
    std::cout << "\n dbar gamma_0 d \n" << std::endl;
    std::cout << C3p[1] << std::endl;
    fout << C3p[1] << std::endl;
    std::cout << "\n ubar O44 u \n" << std::endl;
    std::cout << C3p[2] << std::endl;
    fout << C3p[2] << std::endl;
    std::cout << "\n dbar O44 d \n" << std::endl;
    std::cout << C3p[3] << std::endl;
    fout << C3p[3] << std::endl;

    std::cout << "\n ubar O11 u \n" << std::endl;
    std::cout << C3p[4] << std::endl;
    fout << C3p[4] << std::endl;
    std::cout << "\n dbar O11 d \n" << std::endl;
    std::cout << C3p[5] << std::endl;
    fout << C3p[5] << std::endl;
    std::cout << "\n ubar O22 u \n" << std::endl;
    std::cout << C3p[6] << std::endl;
    fout << C3p[6] << std::endl;
    std::cout << "\n dbar O22 d \n" << std::endl;
    std::cout << C3p[7] << std::endl;
    fout << C3p[7] << std::endl;
    std::cout << "\n ubar O33 u \n" << std::endl;
    std::cout << C3p[8] << std::endl;
    fout << C3p[8] << std::endl;
    std::cout << "\n dbar O33 d \n" << std::endl;
    std::cout << C3p[9] << std::endl;
    fout << C3p[9] << std::endl;

    fout.close();
  }


#ifdef __CALCULATE_TWOPOINT__
  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = size_t(floats["APE_steps"]);
  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = size_t(floats["Jac_steps"]);

  if (weave.isRoot())
  {
    std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
    std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
  }

  Smear::APE APE_tool(APE_alpha);
  Smear::Jacobi Jacobi_tool(Jac_alpha);

  APE_tool.smear(gauge_field, APE_iterations, timeslice_sink);

  forwardProp_u.smearJacobi(Jac_alpha, Jac_iterations, gauge_field, timeslice_sink);
  forwardProp_d.smearJacobi(Jac_alpha, Jac_iterations, gauge_field, timeslice_sink);

  if (weave.isRoot())
    std::cout << "propagators and gauge field smeared successfully\n" << std::endl;

  Core::Correlator C2_P = Contract::proton_twopoint(forwardProp_u, forwardProp_u, forwardProp_d, Base::proj_PARITY_PLUS_TM);


  if (weave.isRoot())
  {
    std::ofstream fout("output_2point.dat");
    std::cout << "proton two point" << std::endl;
    std::cout << C2_P << std::endl;
  }


#endif

  if (weave.isRoot())
    std::cout << "contractions performed and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
