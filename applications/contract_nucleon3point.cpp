// $Id$




// ################################################################################
// ## important NOTE:: we don't want to have all operators with all projectors ####
// ## => use one std::vector< Base::Operator > for each projector #################
// ################################################################################



/*
  This code loads two forward propagators (u,d) and 8 backward propagators (4 different projectors for u and d)
  and calculates nucleon threepoint functions.
  Some modifications might have to be done here concerning which correlation function is calculated and
   what the output is supposed to look like.
*/

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
// for forward and backward propagators separately
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__

// comment this if you don't want the threepoint function to be calculated
#define __CALCULATE_THREEPOINT__

// comment this if you don't want the twopoint function to be calculated
#define __CALCULATE_TWOPOINT__


int main(int argc, char **argv)
{

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./contract_nucleon3point_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::map< std::string, int > operators;
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
  std::vector< std::string > const &seqPropFilesD_proj0(files[ 3]);
  std::vector< std::string > const &seqPropFilesU_proj0(files[ 4]);
  std::vector< std::string > const &seqPropFilesD_proj1(files[ 5]);
  std::vector< std::string > const &seqPropFilesU_proj1(files[ 6]);
  std::vector< std::string > const &seqPropFilesD_proj2(files[ 7]);
  std::vector< std::string > const &seqPropFilesU_proj2(files[ 8]);
  std::vector< std::string > const &seqPropFilesD_proj3(files[ 9]);
  std::vector< std::string > const &seqPropFilesU_proj3(files[10]);

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

  Core::Propagator backwardProp_u_proj0(L, T);
  Core::Propagator backwardProp_d_proj0(L, T);

  Tool::IO::load(&backwardProp_u_proj0, seqPropFilesU_proj0, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark backward propagator (projector 0) successfully loaded\n" << std::endl;

  Tool::IO::load(&backwardProp_d_proj0, seqPropFilesD_proj0, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark backward propagator (projector 0)  successfully loaded\n" << std::endl;


  Core::Propagator backwardProp_u_proj1(L, T);
  Core::Propagator backwardProp_d_proj1(L, T);

  Tool::IO::load(&backwardProp_u_proj1, seqPropFilesU_proj1, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark backward propagator (projector 1) successfully loaded\n" << std::endl;

  Tool::IO::load(&backwardProp_d_proj1, seqPropFilesD_proj1, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark backward propagator (projector 1)  successfully loaded\n" << std::endl;

  Core::Propagator backwardProp_u_proj2(L, T);
  Core::Propagator backwardProp_d_proj2(L, T);

  Tool::IO::load(&backwardProp_u_proj2, seqPropFilesU_proj2, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark backward propagator (projector 2) successfully loaded\n" << std::endl;

  Tool::IO::load(&backwardProp_d_proj2, seqPropFilesD_proj2, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark backward propagator (projector 2)  successfully loaded\n" << std::endl;

  Core::Propagator backwardProp_u_proj3(L, T);
  Core::Propagator backwardProp_d_proj3(L, T);

  Tool::IO::load(&backwardProp_u_proj3, seqPropFilesU_proj3, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark backward propagator (projector 3) successfully loaded\n" << std::endl;

  Tool::IO::load(&backwardProp_d_proj3, seqPropFilesD_proj3, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark backward propagator (projector 3)  successfully loaded\n" << std::endl;



#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__
  backwardProp_u_proj0.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj0.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_u_proj1.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj1.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_u_proj2.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj2.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_u_proj3.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj3.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
#endif

//   // twopoint test
//   {
//     std::vector< Base::Operator > operators;
// 
//     operators.push_back(Base::op_UNITY);
// 
//     std::vector< Core::BaryonCorrelator > C3p =
//       Contract::proton_threepoint_sequential(backwardProp_u_proj0, forwardProp_u,
//                                              backwardProp_d_proj0, forwardProp_d,
//                                              &gauge_field,
//                                              operators);
// 
// 
//   }

  std::vector< Base::Operator > my_operators;

  my_operators.push_back(Base::op_GAMMA_1);
  my_operators.push_back(Base::op_GAMMA_2);
  my_operators.push_back(Base::op_GAMMA_3);
  my_operators.push_back(Base::op_GAMMA_4);
  my_operators.push_back(Base::op_GAMMA_15);
  my_operators.push_back(Base::op_GAMMA_25);
  my_operators.push_back(Base::op_GAMMA_35);
  my_operators.push_back(Base::op_GAMMA_45);
  my_operators.push_back(Base::op_O11);
  my_operators.push_back(Base::op_O22);
  my_operators.push_back(Base::op_O33);
  my_operators.push_back(Base::op_O44);


  //my_operators.push_back(Base::op_CONSERVED_GAMMA_4);

  if (weave.isRoot())
    std::cout << "\n calculating 3-point function(s) \n" << std::endl;


// NOTE:: we don't want to have all operators with all projectors
// => use one std::vector< Base::Operator > for each projector


  std::vector< Core::BaryonCorrelator > C3p_proj0
    = Contract::proton_threepoint_sequential(backwardProp_u_proj0, forwardProp_u,
                                             backwardProp_d_proj0, forwardProp_d,
                                             &gauge_field, my_operators);

  std::vector< Core::BaryonCorrelator > C3p_proj1
    = Contract::proton_threepoint_sequential(backwardProp_u_proj1, forwardProp_u,
                                             backwardProp_d_proj1, forwardProp_d,
                                             &gauge_field, my_operators);

  std::vector< Core::BaryonCorrelator > C3p_proj2
    = Contract::proton_threepoint_sequential(backwardProp_u_proj2, forwardProp_u,
                                             backwardProp_d_proj2, forwardProp_d,
                                             &gauge_field, my_operators);

  std::vector< Core::BaryonCorrelator > C3p_proj3
    = Contract::proton_threepoint_sequential(backwardProp_u_proj3, forwardProp_u,
                                             backwardProp_d_proj3, forwardProp_d,
                                             &gauge_field, my_operators);


  for(size_t i=0; i<2*my_operators.size(); i++)
  {
    (C3p_proj0[i]).setOffset(timeslice_source);
    (C3p_proj1[i]).setOffset(timeslice_source);
    (C3p_proj2[i]).setOffset(timeslice_source);
    (C3p_proj3[i]).setOffset(timeslice_source);
  }


  if (weave.isRoot())
  {
    std::ofstream fout;

    std::vector< std::string > output_filenames_uu;
    std::vector< std::string > output_filenames_dd;

    output_filenames_uu.push_back("output_3point_vector_0_uu.dat");
    output_filenames_uu.push_back("output_3point_axial_0_uu.dat");
    output_filenames_uu.push_back("output_3point_derivative_0_uu.dat");

    output_filenames_dd.push_back("output_3point_vector_0_dd.dat");
    output_filenames_dd.push_back("output_3point_axial_0_dd.dat");
    output_filenames_dd.push_back("output_3point_derivative_0_dd.dat");

    output_filenames_uu.push_back("output_3point_vector_1_uu.dat");
    output_filenames_uu.push_back("output_3point_axial_1_uu.dat");
    output_filenames_uu.push_back("output_3point_derivative_1_uu.dat");

    output_filenames_dd.push_back("output_3point_vector_1_dd.dat");
    output_filenames_dd.push_back("output_3point_axial_1_dd.dat");
    output_filenames_dd.push_back("output_3point_derivative_1_dd.dat");

    output_filenames_uu.push_back("output_3point_vector_2_uu.dat");
    output_filenames_uu.push_back("output_3point_axial_2_uu.dat");
    output_filenames_uu.push_back("output_3point_derivative_2_uu.dat");

    output_filenames_dd.push_back("output_3point_vector_2_dd.dat");
    output_filenames_dd.push_back("output_3point_axial_2_dd.dat");
    output_filenames_dd.push_back("output_3point_derivative_2_dd.dat");

    output_filenames_uu.push_back("output_3point_vector_3_uu.dat");
    output_filenames_uu.push_back("output_3point_axial_3_uu.dat");
    output_filenames_uu.push_back("output_3point_derivative_3_uu.dat");

    output_filenames_dd.push_back("output_3point_vector_3_dd.dat");
    output_filenames_dd.push_back("output_3point_axial_3_dd.dat");
    output_filenames_dd.push_back("output_3point_derivative_3_dd.dat");

    // projector 0

    fout.open(output_filenames_uu[0].c_str());
    fout << C3p_proj0[0] << std::endl;
    fout << C3p_proj0[2] << std::endl;
    fout << C3p_proj0[4] << std::endl;
    fout << C3p_proj0[6] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[0].c_str());
    fout << C3p_proj0[1] << std::endl;
    fout << C3p_proj0[3] << std::endl;
    fout << C3p_proj0[5] << std::endl;
    fout << C3p_proj0[7] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[1].c_str());
    fout << C3p_proj0[ 8] << std::endl;
    fout << C3p_proj0[10] << std::endl;
    fout << C3p_proj0[12] << std::endl;
    fout << C3p_proj0[14] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[1].c_str());
    fout << C3p_proj0[ 9] << std::endl;
    fout << C3p_proj0[11] << std::endl;
    fout << C3p_proj0[13] << std::endl;
    fout << C3p_proj0[15] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[2].c_str());
    fout << C3p_proj0[16] << std::endl;
    fout << C3p_proj0[18] << std::endl;
    fout << C3p_proj0[20] << std::endl;
    fout << C3p_proj0[22] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[2].c_str());
    fout << C3p_proj0[17] << std::endl;
    fout << C3p_proj0[19] << std::endl;
    fout << C3p_proj0[21] << std::endl;
    fout << C3p_proj0[23] << std::endl;
    fout.close();


    // projector 1

    fout.open(output_filenames_uu[3].c_str());
    fout << C3p_proj1[0] << std::endl;
    fout << C3p_proj1[2] << std::endl;
    fout << C3p_proj1[4] << std::endl;
    fout << C3p_proj1[6] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[3].c_str());
    fout << C3p_proj1[1] << std::endl;
    fout << C3p_proj1[3] << std::endl;
    fout << C3p_proj1[5] << std::endl;
    fout << C3p_proj1[7] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[4].c_str());
    fout << C3p_proj1[ 8] << std::endl;
    fout << C3p_proj1[10] << std::endl;
    fout << C3p_proj1[12] << std::endl;
    fout << C3p_proj1[14] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[4].c_str());
    fout << C3p_proj1[ 9] << std::endl;
    fout << C3p_proj1[11] << std::endl;
    fout << C3p_proj1[13] << std::endl;
    fout << C3p_proj1[15] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[5].c_str());
    fout << C3p_proj1[16] << std::endl;
    fout << C3p_proj1[18] << std::endl;
    fout << C3p_proj1[20] << std::endl;
    fout << C3p_proj1[22] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[5].c_str());
    fout << C3p_proj1[17] << std::endl;
    fout << C3p_proj1[19] << std::endl;
    fout << C3p_proj1[21] << std::endl;
    fout << C3p_proj1[23] << std::endl;
    fout.close();


    // projector 2

    fout.open(output_filenames_uu[6].c_str());
    fout << C3p_proj2[0] << std::endl;
    fout << C3p_proj2[2] << std::endl;
    fout << C3p_proj2[4] << std::endl;
    fout << C3p_proj2[6] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[6].c_str());
    fout << C3p_proj2[1] << std::endl;
    fout << C3p_proj2[3] << std::endl;
    fout << C3p_proj2[5] << std::endl;
    fout << C3p_proj2[7] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[7].c_str());
    fout << C3p_proj2[ 8] << std::endl;
    fout << C3p_proj2[10] << std::endl;
    fout << C3p_proj2[12] << std::endl;
    fout << C3p_proj2[14] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[7].c_str());
    fout << C3p_proj2[ 9] << std::endl;
    fout << C3p_proj2[11] << std::endl;
    fout << C3p_proj2[13] << std::endl;
    fout << C3p_proj2[15] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[8].c_str());
    fout << C3p_proj2[16] << std::endl;
    fout << C3p_proj2[18] << std::endl;
    fout << C3p_proj2[20] << std::endl;
    fout << C3p_proj2[22] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[8].c_str());
    fout << C3p_proj2[17] << std::endl;
    fout << C3p_proj2[19] << std::endl;
    fout << C3p_proj2[21] << std::endl;
    fout << C3p_proj2[23] << std::endl;
    fout.close();


    // projector 3

    fout.open(output_filenames_uu[9].c_str());
    fout << C3p_proj3[0] << std::endl;
    fout << C3p_proj3[2] << std::endl;
    fout << C3p_proj3[4] << std::endl;
    fout << C3p_proj3[6] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[9].c_str());
    fout << C3p_proj3[1] << std::endl;
    fout << C3p_proj3[3] << std::endl;
    fout << C3p_proj3[5] << std::endl;
    fout << C3p_proj3[7] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[10].c_str());
    fout << C3p_proj3[ 8] << std::endl;
    fout << C3p_proj3[10] << std::endl;
    fout << C3p_proj3[12] << std::endl;
    fout << C3p_proj3[14] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[10].c_str());
    fout << C3p_proj3[ 9] << std::endl;
    fout << C3p_proj3[11] << std::endl;
    fout << C3p_proj3[13] << std::endl;
    fout << C3p_proj3[15] << std::endl;
    fout.close();

    fout.open(output_filenames_uu[11].c_str());
    fout << C3p_proj3[16] << std::endl;
    fout << C3p_proj3[18] << std::endl;
    fout << C3p_proj3[20] << std::endl;
    fout << C3p_proj3[22] << std::endl;
    fout.close();
    fout.open(output_filenames_dd[11].c_str());
    fout << C3p_proj3[17] << std::endl;
    fout << C3p_proj3[19] << std::endl;
    fout << C3p_proj3[21] << std::endl;
    fout << C3p_proj3[23] << std::endl;
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
