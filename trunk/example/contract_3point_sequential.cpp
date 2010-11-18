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
// for forward and backward propagators separately
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__

// comment this if you don't want the threepoint function to be calculated
#define __CALCULATE_THREEPOINT__

// comment this if you don't want the twopoint function to be calculated
#define __CALCULATE_TWOPOINT__

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

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./contract_3point_sequential_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::vector< int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);

  size_t const L(L_tmp);
  size_t const T(T_tmp);

  Base::Weave weave(L, T);

  if (weave.isRoot())
    std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu    = floats["mu"];
  if (weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

  size_t const sourceSinkSeparation = size_t(floats["sourceSinkSeparation"]);
  assert(sourceSinkSeparation > 0 && sourceSinkSeparation < T-1);


  size_t const * const source_position = positions[0];
  size_t const timeslice_source = source_position[Base::idx_T] % T;
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  size_t const timeslice_sink = (timeslice_source +  sourceSinkSeparation) % T;
  if (weave.isRoot())
    std::cout << "timeslice (sink) = " << timeslice_sink << std::endl;

  // make sure the boundary is not crossed by source-sink correlaton function
  size_t const timeslice_boundary = (timeslice_source + (T/2)) % T;
  if (weave.isRoot())
    std::cout << "timeslice (boundary) = " << timeslice_boundary << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
  std::vector< std::string > const &seqPropFilesD(files[3]);
  std::vector< std::string > const &seqPropFilesU(files[4]);

  Core::Field< QCD::Gauge > gauge_field(L, T);
  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;

  Core::Propagator forwardProp_u(L, T);
  Tool::IO::load(&forwardProp_u, propfilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark forward propagator successfully loaded\n" << std::endl;

  Core::Propagator forwardProp_d(L, T);
  Tool::IO::load(&forwardProp_d, propfilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark forward propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
  forwardProp_d.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  forwardProp_u.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
#endif


#ifdef __CALCULATE_THREEPOINT__

  Core::Propagator backwardProp_u(L, T);
  Core::Propagator backwardProp_d(L, T);

  Tool::IO::load(&backwardProp_u, seqPropFilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark backward propagator successfully loaded\n" << std::endl;

  Tool::IO::load(&backwardProp_d, seqPropFilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark backward propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__
  backwardProp_u.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
#endif


  std::vector< Base::Operator > my_operators;

//   my_operators.push_back(Base::op_GAMMA_4);
//   my_operators.push_back(Base::op_GAMMA_1);
//   my_operators.push_back(Base::op_GAMMA_2);
//   my_operators.push_back(Base::op_GAMMA_3);
#ifdef __GAMMA_1__
  my_operators.push_back(Base::op_GAMMA_15);
#endif
#ifdef __GAMMA_2__
  my_operators.push_back(Base::op_GAMMA_25);
#endif
#ifdef __GAMMA_3__
  my_operators.push_back(Base::op_GAMMA_35);
#endif
#ifdef __GAMMA_4__
  my_operators.push_back(Base::op_GAMMA_45);
#endif
//   my_operators.push_back(Base::op_O44);
//   my_operators.push_back(Base::op_O11);
//   my_operators.push_back(Base::op_O22);
//   my_operators.push_back(Base::op_O33);
//   my_operators.push_back(Base::op_CONSERVED_GAMMA_4);

  if (weave.isRoot())
    std::cout << "\n calculating 3-point function(s) \n" << std::endl;


  std::vector< Core::BaryonCorrelator > C3p = Contract::proton_threepoint_sequential(backwardProp_u, forwardProp_u,
                                                                               backwardProp_d, forwardProp_d,
                                                                               &gauge_field,
                                                                               my_operators);

  for(size_t i=0; i<2*my_operators.size(); i++)
    (C3p[i]).setOffset(timeslice_source);

  if (weave.isRoot())
  {
   #ifdef __GAMMA_1__
   std::ofstream fout("output_3point_axial_1_uu.dat");
   #endif
   #ifdef __GAMMA_2__
   std::ofstream fout("output_3point_axial_2_uu.dat");
   #endif
   #ifdef __GAMMA_3__
   std::ofstream fout("output_3point_axial_3_uu.dat");
   #endif
   #ifdef __GAMMA_4__
   std::ofstream fout("output_3point_axial_4_uu.dat");
   #endif
   fout << C3p[0] << std::endl;
   fout.close();
   #ifdef __GAMMA_1__
   fout.open("output_3point_axial_1_dd.dat");
   #endif
   #ifdef __GAMMA_2__
   fout.open("output_3point_axial_2_dd.dat");
   #endif
   #ifdef __GAMMA_3__
   fout.open("output_3point_axial_3_dd.dat");
   #endif
   #ifdef __GAMMA_4__
   fout.open("output_3point_axial_4_dd.dat");
   #endif
   fout << C3p[1] << std::endl;
   fout.close();
  }

#endif

#ifdef __CALCULATE_TWOPOINT__

  if (weave.isRoot())
    std::cout << "\n calculating 2-point function \n" << std::endl;

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

  APE_tool.smear(gauge_field, APE_iterations);

  forwardProp_u.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  forwardProp_d.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  if (weave.isRoot())
    std::cout << "propagators and gauge field smeared successfully\n" << std::endl;

  Core::BaryonCorrelator C2_P = Contract::proton_twopoint(forwardProp_u, forwardProp_u, forwardProp_d, Base::proj_PARITY_PLUS_TM);
  C2_P.setOffset(timeslice_source);

  if (weave.isRoot())
  {
    std::ofstream fout("output_2point.dat");
    std::cout << "proton two point" << std::endl;
    std::cout << C2_P << std::endl;
    fout << C2_P << std::endl;
    fout.close();
  }


#endif

  if (weave.isRoot())
    std::cout << "contractions performed and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
