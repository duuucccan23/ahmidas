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

// do we want to use all 4 projectors? (otherwise onlu (1+gamma_0)/2 is used)
#define __4_PROJECTORS__


int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);
  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./contract_nucleon3point_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::vector< int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);

  size_t const L(L_tmp);
  size_t const T(T_tmp);

  Base::Weave weave(L, T);

  std::vector< int* > momenta;
  for(size_t I=0; I<257; I++)
    momenta.push_back( new int[3]);
  {
    int momenta_raw[771] = {
    +0, +0, +0,
     1, +0, +0,
    -1, +0, +0,
    +0,  1, +0,
    +0, -1, +0,
    +0, +0,  1,
    +0, +0, -1,
     1, +0,  1,
    -1, +0, -1,
     1, +0, -1,
    -1, +0,  1,
     1,  1, +0,
    -1, -1, +0,
     1, -1, +0,
    -1,  1, +0,
    +0,  1,  1,
    +0, -1, -1,
    +0,  1, -1,
    +0, -1,  1,
     1,  1,  1,
    -1, -1, -1,
     1,  1, -1,
    -1, -1,  1,
     1, -1,  1,
    -1,  1, -1,
     1, -1, -1,
    -1,  1,  1,
     2, +0, +0,
    -2, +0, +0,
    +0,  2, +0,
    +0, -2, +0,
    +0, +0,  2,
    +0, +0, -2,
     2,  1, +0,
    -2, -1, +0,
     2, -1, +0,
    -2,  1, +0,
     2, +0,  1,
    -2, +0, -1,
     2, +0, -1,
    -2, +0,  1,
     1,  2, +0,
    -1, -2, +0,
     1, -2, +0,
    -1,  2, +0,
     1, +0,  2,
    -1, +0, -2,
     1, +0, -2,
    -1, +0,  2,
    +0,  2,  1,
    +0, -2, -1,
    +0,  2, -1,
    +0, -2,  1,
    +0,  1,  2,
    +0, -1, -2,
    +0,  1, -2,
    +0, -1,  2,
     2,  1,  1,
    -2, -1, -1,
     2,  1, -1,
    -2, -1,  1,
     2, -1,  1,
    -2,  1, -1,
     2, -1, -1,
    -2,  1,  1,
     1,  2,  1,
    -1, -2, -1,
     1,  2, -1,
    -1, -2,  1,
    -1,  2,  1,
     1, -2, -1,
    -1,  2, -1,
     1, -2,  1,
     1,  1,  2,
    -1, -1, -2,
     1, -1,  2,
    -1,  1, -2,
    -1,  1,  2,
     1, -1, -2,
    -1, -1,  2,
     1,  1, -2,
     2, +0,  2,
    -2, +0, -2,
     2, +0, -2,
    -2, +0,  2,
     2,  2, +0,
    -2, -2, +0,
     2, -2, +0,
    -2,  2, +0,
    +0,  2,  2,
    +0, -2, -2,
    +0,  2, -2,
    +0, -2,  2,
     3, +0, +0,
    -3, +0, +0,
    +0,  3, +0,
    +0, -3, +0,
    +0, +0,  3,
    +0, +0, -3,
     1,  2,  2,
    -1, -2, -2,
     1,  2, -2,
    -1, -2,  2,
     1, -2,  2,
    -1,  2, -2,
     1, -2, -2,
    -1,  2,  2,
     2,  1,  2,
    -2, -1, -2,
     2,  1, -2,
    -2, -1,  2,
    -2,  1,  2,
     2, -1, -2,
    -2,  1, -2,
     2, -1,  2,
     2,  2,  1,
    -2, -2, -1,
     2, -2,  1,
    -2,  2, -1,
    -2,  2,  1,
     2, -2, -1,
    -2, -2,  1,
     2,  2, -1,
     3,  1, +0,
    -3, -1, +0,
     3, -1, +0,
    -3,  1, +0,
     3, +0,  1,
    -3, +0, -1,
     3, +0, -1,
    -3, +0,  1,
     1,  3, +0,
    -1, -3, +0,
     1, -3, +0,
    -1,  3, +0,
     1, +0,  3,
    -1, +0, -3,
     1, +0, -3,
    -1, +0,  3,
    +0,  3,  1,
    +0, -3, -1,
    +0,  3, -1,
    +0, -3,  1,
    +0,  1,  3,
    +0, -1, -3,
    +0,  1, -3,
    +0, -1,  3,
     3,  1,  1,
    -3, -1, -1,
     3,  1, -1,
    -3, -1,  1,
     3, -1,  1,
    -3,  1, -1,
     3, -1, -1,
    -3,  1,  1,
     1,  3,  1,
    -1, -3, -1,
     1,  3, -1,
    -1, -3,  1,
    -1,  3,  1,
     1, -3, -1,
    -1,  3, -1,
     1, -3,  1,
     1,  1,  3,
    -1, -1, -3,
     1, -1,  3,
    -1,  1, -3,
    -1,  1,  3,
     1, -1, -3,
    -1, -1,  3,
     1,  1, -3,
     2,  2,  2,
    -2, -2, -2,
     2,  2, -2,
    -2, -2,  2,
     2, -2,  2,
    -2,  2, -2,
     2, -2, -2,
    -2,  2,  2,
     3,  2, +0,
    -3, -2, +0,
     3, -2, +0,
    -3,  2, +0,
     3, +0,  2,
    -3, +0, -2,
     3, +0, -2,
    -3, +0,  2,
     2,  3, +0,
    -2, -3, +0,
     2, -3, +0,
    -2,  3, +0,
     2, +0,  3,
    -2, +0, -3,
     2, +0, -3,
    -2, +0,  3,
    +0,  3,  2,
    +0, -3, -2,
    +0,  3, -2,
    +0, -3,  2,
    +0,  2,  3,
    +0, -2, -3,
    +0,  2, -3,
    +0, -2,  3,
     3,  2,  1,
    -3, -2, -1,
    -3,  2,  1,
     3, -2, -1,
     3, -2,  1,
    -3,  2, -1,
     3,  2, -1,
    -3, -2,  1,
     2,  3,  1,
    -2, -3, -1,
    -2,  3,  1,
     2, -3, -1,
     2, -3,  1,
    -2,  3, -1,
     2,  3, -1,
    -2, -3,  1,
     3,  1,  2,
    -3, -1, -2,
    -3,  1,  2,
     3, -1, -2,
     3, -1,  2,
    -3,  1, -2,
     3,  1, -2,
    -3, -1,  2,
     2,  1,  3,
    -2, -1, -3,
    -2,  1,  3,
     2, -1, -3,
     2, -1,  3,
    -2,  1, -3,
     2,  1, -3,
    -2, -1,  3,
     1,  2,  3,
    -1, -2, -3,
    -1,  2,  3,
     1, -2, -3,
     1, -2,  3,
    -1,  2, -3,
     1,  2, -3,
    -1, -2,  3,
     1,  3,  2,
    -1, -3, -2,
    -1,  3,  2,
     1, -3, -2,
     1, -3,  2,
    -1,  3, -2,
     1,  3, -2,
    -1, -3,  2,
     4, +0, +0,
    -4, +0, +0,
    +0,  4, +0,
    +0, -4, +0,
    +0, +0,  4,
    +0, +0, -4};

    for(size_t I=0; I<momenta.size(); I++)
      std::copy(&(momenta_raw[3*I]), &(momenta_raw[3*I]) + 3, momenta[I]);
  }

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
//   size_t const timeslice_boundary = (timeslice_source + (T/2)) % T;
  size_t const timeslice_boundary = (timeslice_sink + 2) % T;
  if (weave.isRoot())
    std::cout << "timeslice (boundary) = " << timeslice_boundary << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
  std::vector< std::string > const &seqPropFilesD_proj0(files[ 3]);
  std::vector< std::string > const &seqPropFilesU_proj0(files[ 4]);

#ifdef __4_PROJECTORS__
  std::vector< std::string > const &seqPropFilesD_proj1(files[ 5]);
  std::vector< std::string > const &seqPropFilesU_proj1(files[ 6]);
  std::vector< std::string > const &seqPropFilesD_proj2(files[ 7]);
  std::vector< std::string > const &seqPropFilesU_proj2(files[ 8]);
  std::vector< std::string > const &seqPropFilesD_proj3(files[ 9]);
  std::vector< std::string > const &seqPropFilesU_proj3(files[10]);
#endif

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

#ifdef __4_PROJECTORS__

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
#endif


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__
  backwardProp_u_proj0.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj0.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
 #ifdef __4_PROJECTORS__
  backwardProp_u_proj1.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj1.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_u_proj2.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj2.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_u_proj3.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  backwardProp_d_proj3.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  #endif
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

  //my_operators.push_back(Base::op_CONSERVED_GAMMA_4);

  if (weave.isRoot())
    std::cout << "\n calculating 3-point function(s) \n" << std::endl;


  if (weave.isRoot())
    std::cout << "\n local currents ... \n" << std::endl;


  // ********************************************************************************************
  // local currents *****************************************************************************
  // ********************************************************************************************
  {
  std::vector< Core::BaryonCorrelator > C3p_proj0_local
    = Contract::proton_threepoint_sequential(backwardProp_u_proj0, forwardProp_u,
                                             backwardProp_d_proj0, forwardProp_d,
                                             &gauge_field, "local");

  #ifdef __4_PROJECTORS__

  std::vector< Core::BaryonCorrelator > C3p_proj1_local
    = Contract::proton_threepoint_sequential(backwardProp_u_proj1, forwardProp_u,
                                             backwardProp_d_proj1, forwardProp_d,
                                             &gauge_field, "local");
  std::vector< Core::BaryonCorrelator > C3p_proj2_local
    = Contract::proton_threepoint_sequential(backwardProp_u_proj2, forwardProp_u,
                                             backwardProp_d_proj2, forwardProp_d,
                                             &gauge_field, "local");
  std::vector< Core::BaryonCorrelator > C3p_proj3_local
      = Contract::proton_threepoint_sequential(backwardProp_u_proj3, forwardProp_u,
                                               backwardProp_d_proj3, forwardProp_d,
                                               &gauge_field, "local");

  #endif


  int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
  C3p_proj0_local[0].prepareMomentumProjection(sourcePos);

  std::ofstream *fout_proj0_uu = NULL;
  std::ofstream *fout_proj0_dd = NULL;
  #ifdef __4_PROJECTORS__
  std::ofstream *fout_proj1_uu = NULL;
  std::ofstream *fout_proj1_dd = NULL;
  std::ofstream *fout_proj2_uu = NULL;
  std::ofstream *fout_proj2_dd = NULL;
  std::ofstream *fout_proj3_uu = NULL;
  std::ofstream *fout_proj3_dd = NULL;
  #endif


  if (weave.isRoot())
  {
    fout_proj0_uu =  new std::ofstream("output_3point_local_proj0_uu.dat");
    fout_proj0_dd =  new std::ofstream("output_3point_local_proj0_dd.dat");
    #ifdef __4_PROJECTORS__
    fout_proj1_uu =  new std::ofstream("output_3point_local_proj1_uu.dat");
    fout_proj1_dd =  new std::ofstream("output_3point_local_proj1_dd.dat");
    fout_proj2_uu =  new std::ofstream("output_3point_local_proj2_uu.dat");
    fout_proj2_dd =  new std::ofstream("output_3point_local_proj2_dd.dat");
    fout_proj3_uu =  new std::ofstream("output_3point_local_proj3_uu.dat");
    fout_proj3_dd =  new std::ofstream("output_3point_local_proj3_dd.dat");
    #endif
  }


  assert(C3p_proj0_local.size() == 32);

  for(size_t i=0; i<C3p_proj0_local.size(); i++)
  {
    std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_local[i].momentumProjection(momenta));
    #ifdef __4_PROJECTORS__
    std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_local[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_local[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_local[i].momentumProjection(momenta));
    #endif


    for(size_t I=0; I<momenta.size(); I++)
    {
      std::ostringstream oss;
      oss << std::setw(3) << i/2 << " " << std::flush;
      std::string const prefix(oss.str());
      all_corrs_proj0[I].setOffset(timeslice_source);
      #ifdef __4_PROJECTORS__
      all_corrs_proj1[I].setOffset(timeslice_source);
      all_corrs_proj2[I].setOffset(timeslice_source);
      all_corrs_proj3[I].setOffset(timeslice_source);
      #endif
      if (weave.isRoot())
      {
        if(i%2 == 0)
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          #endif
        }
        else
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          #endif
        }
      }
    }

    weave.barrier();

  }

  if (weave.isRoot())
  {
    fout_proj0_uu->close();
    fout_proj0_dd->close();
    #ifdef __4_PROJECTORS__
    fout_proj1_uu->close();
    fout_proj1_dd->close();
    fout_proj2_uu->close();
    fout_proj2_dd->close();
    fout_proj3_uu->close();
    fout_proj3_dd->close();
    #endif
  }

  delete fout_proj0_uu;
  delete fout_proj0_dd;
  #ifdef __4_PROJECTORS__
  delete fout_proj1_uu;
  delete fout_proj1_dd;
  delete fout_proj2_uu;
  delete fout_proj2_dd;
  delete fout_proj3_uu;
  delete fout_proj3_dd;
  #endif
  }


  if (weave.isRoot())
    std::cout << "\n conserved currents ... \n" << std::endl;


  // ********************************************************************************************
  // conserved currents *************************************************************************
  // ********************************************************************************************
  {
  std::vector< Core::BaryonCorrelator > C3p_proj0_conserved_axial;
  std::vector< Core::BaryonCorrelator > C3p_proj0_conserved_vector
    = Contract::proton_threepoint_sequential(backwardProp_u_proj0, forwardProp_u,
                                             backwardProp_d_proj0, forwardProp_d,
                                             &gauge_field, "noether", &C3p_proj0_conserved_axial);

  #ifdef __4_PROJECTORS__
  std::vector< Core::BaryonCorrelator > C3p_proj1_conserved_axial;
  std::vector< Core::BaryonCorrelator > C3p_proj1_conserved_vector
    = Contract::proton_threepoint_sequential(backwardProp_u_proj1, forwardProp_u,
                                             backwardProp_d_proj1, forwardProp_d,
                                             &gauge_field, "noether", &C3p_proj1_conserved_axial);
  std::vector< Core::BaryonCorrelator > C3p_proj2_conserved_axial;
  std::vector< Core::BaryonCorrelator > C3p_proj2_conserved_vector
    = Contract::proton_threepoint_sequential(backwardProp_u_proj2, forwardProp_u,
                                             backwardProp_d_proj2, forwardProp_d,
                                             &gauge_field, "noether", &C3p_proj2_conserved_axial);
  std::vector< Core::BaryonCorrelator > C3p_proj3_conserved_axial;
  std::vector< Core::BaryonCorrelator > C3p_proj3_conserved_vector
    = Contract::proton_threepoint_sequential(backwardProp_u_proj3, forwardProp_u,
                                             backwardProp_d_proj3, forwardProp_d,
                                             &gauge_field, "noether", &C3p_proj3_conserved_axial);

  #endif


  int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
  C3p_proj0_conserved_vector[0].prepareMomentumProjection(sourcePos);

  std::ofstream *fout_proj0_uu = NULL;
  std::ofstream *fout_proj0_dd = NULL;
  #ifdef __4_PROJECTORS__
  std::ofstream *fout_proj1_uu = NULL;
  std::ofstream *fout_proj1_dd = NULL;
  std::ofstream *fout_proj2_uu = NULL;
  std::ofstream *fout_proj2_dd = NULL;
  std::ofstream *fout_proj3_uu = NULL;
  std::ofstream *fout_proj3_dd = NULL;
  #endif


  if (weave.isRoot())
  {
    fout_proj0_uu =  new std::ofstream("output_3point_conserved_vector_proj0_uu.dat");
    fout_proj0_dd =  new std::ofstream("output_3point_conserved_vector_proj0_dd.dat");
    #ifdef __4_PROJECTORS__
    fout_proj1_uu =  new std::ofstream("output_3point_conserved_vector_proj1_uu.dat");
    fout_proj1_dd =  new std::ofstream("output_3point_conserved_vector_proj1_dd.dat");
    fout_proj2_uu =  new std::ofstream("output_3point_conserved_vector_proj2_uu.dat");
    fout_proj2_dd =  new std::ofstream("output_3point_conserved_vector_proj2_dd.dat");
    fout_proj3_uu =  new std::ofstream("output_3point_conserved_vector_proj3_uu.dat");
    fout_proj3_dd =  new std::ofstream("output_3point_conserved_vector_proj3_dd.dat");
    #endif
  }


  assert(C3p_proj0_conserved_vector.size() == 8);

  for(size_t i=0; i<C3p_proj0_conserved_vector.size(); i++)
  {
    std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_conserved_vector[i].momentumProjection(momenta));
    #ifdef __4_PROJECTORS__
    std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_conserved_vector[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_conserved_vector[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_conserved_vector[i].momentumProjection(momenta));
    #endif


    for(size_t I=0; I<momenta.size(); I++)
    {
      std::ostringstream oss;
      oss << std::setw(3) << i/2 << " " << std::flush;
      std::string const prefix(oss.str());
      all_corrs_proj0[I].setOffset(timeslice_source);
      #ifdef __4_PROJECTORS__
      all_corrs_proj1[I].setOffset(timeslice_source);
      all_corrs_proj2[I].setOffset(timeslice_source);
      all_corrs_proj3[I].setOffset(timeslice_source);
      #endif
      if (weave.isRoot())
      {
        if(i%2 == 0)
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          #endif
        }
        else
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          #endif
        }
      }
    }

    weave.barrier();

  }

  if (weave.isRoot())
  {
    fout_proj0_uu->close();
    fout_proj0_dd->close();
    #ifdef __4_PROJECTORS__
    fout_proj1_uu->close();
    fout_proj1_dd->close();
    fout_proj2_uu->close();
    fout_proj2_dd->close();
    fout_proj3_uu->close();
    fout_proj3_dd->close();
    #endif
  }

  delete fout_proj0_uu;
  delete fout_proj0_dd;
  #ifdef __4_PROJECTORS__
  delete fout_proj1_uu;
  delete fout_proj1_dd;
  delete fout_proj2_uu;
  delete fout_proj2_dd;
  delete fout_proj3_uu;
  delete fout_proj3_dd;
  #endif

  if (weave.isRoot())
  {
    fout_proj0_uu =  new std::ofstream("output_3point_conserved_axial_proj0_uu.dat");
    fout_proj0_dd =  new std::ofstream("output_3point_conserved_axial_proj0_dd.dat");
    #ifdef __4_PROJECTORS__
    fout_proj1_uu =  new std::ofstream("output_3point_conserved_axial_proj1_uu.dat");
    fout_proj1_dd =  new std::ofstream("output_3point_conserved_axial_proj1_dd.dat");
    fout_proj2_uu =  new std::ofstream("output_3point_conserved_axial_proj2_uu.dat");
    fout_proj2_dd =  new std::ofstream("output_3point_conserved_axial_proj2_dd.dat");
    fout_proj3_uu =  new std::ofstream("output_3point_conserved_axial_proj3_uu.dat");
    fout_proj3_dd =  new std::ofstream("output_3point_conserved_axial_proj3_dd.dat");
    #endif
  }


  assert(C3p_proj0_conserved_axial.size() == 8);

  for(size_t i=0; i<C3p_proj0_conserved_axial.size(); i++)
  {
    std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_conserved_axial[i].momentumProjection(momenta));
    #ifdef __4_PROJECTORS__
    std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_conserved_axial[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_conserved_axial[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_conserved_axial[i].momentumProjection(momenta));
    #endif


    for(size_t I=0; I<momenta.size(); I++)
    {
      std::ostringstream oss;
      oss << std::setw(3) << i/2 << " " << std::flush;
      std::string const prefix(oss.str());
      all_corrs_proj0[I].setOffset(timeslice_source);
      #ifdef __4_PROJECTORS__
      all_corrs_proj1[I].setOffset(timeslice_source);
      all_corrs_proj2[I].setOffset(timeslice_source);
      all_corrs_proj3[I].setOffset(timeslice_source);
      #endif
      if (weave.isRoot())
      {
        if(i%2 == 0)
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          #endif
        }
        else
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          #endif
        }
      }
    }
    weave.barrier();
  }

  if (weave.isRoot())
  {
    fout_proj0_uu->close();
    fout_proj0_dd->close();
    #ifdef __4_PROJECTORS__
    fout_proj1_uu->close();
    fout_proj1_dd->close();
    fout_proj2_uu->close();
    fout_proj2_dd->close();
    fout_proj3_uu->close();
    fout_proj3_dd->close();
    #endif
  }

  delete fout_proj0_uu;
  delete fout_proj0_dd;
  #ifdef __4_PROJECTORS__
  delete fout_proj1_uu;
  delete fout_proj1_dd;
  delete fout_proj2_uu;
  delete fout_proj2_dd;
  delete fout_proj3_uu;
  delete fout_proj3_dd;
  #endif
  }


  if (weave.isRoot())
    std::cout << "\n one-derivatives ... \n" << std::endl;

  // ********************************************************************************************
  // one-derivative operators *******************************************************************
  // ********************************************************************************************
  {
  std::vector< Core::BaryonCorrelator > C3p_proj0_1D_polarized;
  std::vector< Core::BaryonCorrelator > C3p_proj0_1D_unpolarized
    = Contract::proton_threepoint_sequential(backwardProp_u_proj0, forwardProp_u,
                                             backwardProp_d_proj0, forwardProp_d,
                                             &gauge_field, "1D", &C3p_proj0_1D_polarized);

  #ifdef __4_PROJECTORS__
  std::vector< Core::BaryonCorrelator > C3p_proj1_1D_polarized;
  std::vector< Core::BaryonCorrelator > C3p_proj1_1D_unpolarized
    = Contract::proton_threepoint_sequential(backwardProp_u_proj1, forwardProp_u,
                                             backwardProp_d_proj1, forwardProp_d,
                                             &gauge_field, "1D", &C3p_proj1_1D_polarized);
  std::vector< Core::BaryonCorrelator > C3p_proj2_1D_polarized;
  std::vector< Core::BaryonCorrelator > C3p_proj2_1D_unpolarized
    = Contract::proton_threepoint_sequential(backwardProp_u_proj2, forwardProp_u,
                                             backwardProp_d_proj2, forwardProp_d,
                                             &gauge_field, "1D", &C3p_proj2_1D_polarized);
  std::vector< Core::BaryonCorrelator > C3p_proj3_1D_polarized;
  std::vector< Core::BaryonCorrelator > C3p_proj3_1D_unpolarized
    = Contract::proton_threepoint_sequential(backwardProp_u_proj3, forwardProp_u,
                                             backwardProp_d_proj3, forwardProp_d,
                                             &gauge_field, "1D", &C3p_proj3_1D_polarized);

  #endif


  int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
  C3p_proj0_1D_unpolarized[0].prepareMomentumProjection(sourcePos);

  std::ofstream *fout_proj0_uu = NULL;
  std::ofstream *fout_proj0_dd = NULL;
  #ifdef __4_PROJECTORS__
  std::ofstream *fout_proj1_uu = NULL;
  std::ofstream *fout_proj1_dd = NULL;
  std::ofstream *fout_proj2_uu = NULL;
  std::ofstream *fout_proj2_dd = NULL;
  std::ofstream *fout_proj3_uu = NULL;
  std::ofstream *fout_proj3_dd = NULL;
  #endif


  if (weave.isRoot())
  {
    fout_proj0_uu =  new std::ofstream("output_3point_1D_unpolarized_proj0_uu.dat");
    fout_proj0_dd =  new std::ofstream("output_3point_1D_unpolarized_proj0_dd.dat");
    #ifdef __4_PROJECTORS__
    fout_proj1_uu =  new std::ofstream("output_3point_1D_unpolarized_proj1_uu.dat");
    fout_proj1_dd =  new std::ofstream("output_3point_1D_unpolarized_proj1_dd.dat");
    fout_proj2_uu =  new std::ofstream("output_3point_1D_unpolarized_proj2_uu.dat");
    fout_proj2_dd =  new std::ofstream("output_3point_1D_unpolarized_proj2_dd.dat");
    fout_proj3_uu =  new std::ofstream("output_3point_1D_unpolarized_proj3_uu.dat");
    fout_proj3_dd =  new std::ofstream("output_3point_1D_unpolarized_proj3_dd.dat");
    #endif
  }


  assert(C3p_proj0_1D_unpolarized.size() == 32);

  for(size_t i=0; i<C3p_proj0_1D_unpolarized.size(); i++)
  {
    std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_1D_unpolarized[i].momentumProjection(momenta));
    #ifdef __4_PROJECTORS__
    std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_1D_unpolarized[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_1D_unpolarized[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_1D_unpolarized[i].momentumProjection(momenta));
    #endif


    for(size_t I=0; I<momenta.size(); I++)
    {
      std::ostringstream oss;
      oss << std::setw(3) << i/2 << " " << std::flush;
      std::string const prefix(oss.str());
      all_corrs_proj0[I].setOffset(timeslice_source);
      #ifdef __4_PROJECTORS__
      all_corrs_proj1[I].setOffset(timeslice_source);
      all_corrs_proj2[I].setOffset(timeslice_source);
      all_corrs_proj3[I].setOffset(timeslice_source);
      #endif
      if (weave.isRoot())
      {
        if(i%2 == 0)
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          #endif
        }
        else
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          #endif
        }
      }
    }

    weave.barrier();

  }

  if (weave.isRoot())
  {
    fout_proj0_uu->close();
    fout_proj0_dd->close();
    #ifdef __4_PROJECTORS__
    fout_proj1_uu->close();
    fout_proj1_dd->close();
    fout_proj2_uu->close();
    fout_proj2_dd->close();
    fout_proj3_uu->close();
    fout_proj3_dd->close();
    #endif
  }

  delete fout_proj0_uu;
  delete fout_proj0_dd;
  #ifdef __4_PROJECTORS__
  delete fout_proj1_uu;
  delete fout_proj1_dd;
  delete fout_proj2_uu;
  delete fout_proj2_dd;
  delete fout_proj3_uu;
  delete fout_proj3_dd;
  #endif

  if (weave.isRoot())
  {
    fout_proj0_uu =  new std::ofstream("output_3point_1D_polarized_proj0_uu.dat");
    fout_proj0_dd =  new std::ofstream("output_3point_1D_polarized_proj0_dd.dat");
    #ifdef __4_PROJECTORS__
    fout_proj1_uu =  new std::ofstream("output_3point_1D_polarized_proj1_uu.dat");
    fout_proj1_dd =  new std::ofstream("output_3point_1D_polarized_proj1_dd.dat");
    fout_proj2_uu =  new std::ofstream("output_3point_1D_polarized_proj2_uu.dat");
    fout_proj2_dd =  new std::ofstream("output_3point_1D_polarized_proj2_dd.dat");
    fout_proj3_uu =  new std::ofstream("output_3point_1D_polarized_proj3_uu.dat");
    fout_proj3_dd =  new std::ofstream("output_3point_1D_polarized_proj3_dd.dat");
    #endif
  }


  assert(C3p_proj0_1D_polarized.size() == 32);

  for(size_t i=0; i<C3p_proj0_1D_polarized.size(); i++)
  {
    std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_1D_polarized[i].momentumProjection(momenta));
    #ifdef __4_PROJECTORS__
    std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_1D_polarized[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_1D_polarized[i].momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_1D_polarized[i].momentumProjection(momenta));
    #endif


    for(size_t I=0; I<momenta.size(); I++)
    {
      std::ostringstream oss;
      oss << std::setw(3) << i/2 << " " << std::flush;
      std::string const prefix(oss.str());
      all_corrs_proj0[I].setOffset(timeslice_source);
      #ifdef __4_PROJECTORS__
      all_corrs_proj1[I].setOffset(timeslice_source);
      all_corrs_proj2[I].setOffset(timeslice_source);
      all_corrs_proj3[I].setOffset(timeslice_source);
      #endif
      if (weave.isRoot())
      {
        if(i%2 == 0)
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          #endif
        }
        else
        {
          all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
          #ifdef __4_PROJECTORS__
          all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
          all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
          all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          #endif
        }
      }
    }
    weave.barrier();
  }

  if (weave.isRoot())
  {
    fout_proj0_uu->close();
    fout_proj0_dd->close();
    #ifdef __4_PROJECTORS__
    fout_proj1_uu->close();
    fout_proj1_dd->close();
    fout_proj2_uu->close();
    fout_proj2_dd->close();
    fout_proj3_uu->close();
    fout_proj3_dd->close();
    #endif
  }

  delete fout_proj0_uu;
  delete fout_proj0_dd;
  #ifdef __4_PROJECTORS__
  delete fout_proj1_uu;
  delete fout_proj1_dd;
  delete fout_proj2_uu;
  delete fout_proj2_dd;
  delete fout_proj3_uu;
  delete fout_proj3_dd;
  #endif
  }

  //***************************************************************************************

#endif
// (#ifdef __CALCULATE_THREEPOINT__)

#ifdef __CALCULATE_TWOPOINT__

  if (weave.isRoot())
    std::cout << "\n calculating 2-point function\n" << std::endl;

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

//   Core::BaryonCorrelator C2_P = Contract::proton_twopoint(forwardProp_u, forwardProp_u, forwardProp_d,
//                                                           Base::proj_PARITY_PLUS_TM);
//   C2_P.setOffset(timeslice_source);
// 
//   if (weave.isRoot())
//   {
//     std::ofstream fout("output_2point.dat");
//     std::cout << "proton two point" << std::endl;
//     std::cout << C2_P << std::endl;
//     fout << C2_P << std::endl;
//     fout.close();
//   }

  {
    if (weave.isRoot())
      std::cout << "\n calculating 2-point function (unprojected) \n" << std::endl;

    forwardProp_u.rotateToPhysicalBasis(true);
    forwardProp_d.rotateToPhysicalBasis(false);

    Core::BaryonCorrelator C2_P = Contract::proton_twopoint(forwardProp_u, forwardProp_u, forwardProp_d, Base::proj_NO_PROJECTOR);
    Core::BaryonCorrelator C2_N = Contract::proton_twopoint(forwardProp_d, forwardProp_d, forwardProp_u, Base::proj_NO_PROJECTOR);

    // C2_P.deleteField();


    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C2_P.prepareMomentumProjection(sourcePos);


    std::vector< Core::BaryonCorrelator > all_corrsP(C2_P.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C2_N.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_2point_proton.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        // all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        // all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    for(size_t I=0; I<momenta.size(); I++)
      delete [] momenta[I];
    momenta.clear();
    delete fout;

    weave.barrier();

  }

#endif
// (#ifdef __CALCULATE_TWOPOINT__)

  if (weave.isRoot())
    std::cout << "contractions performed and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
