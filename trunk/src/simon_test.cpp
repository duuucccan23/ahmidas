#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Ahmidas.h>
#include <L0/Dirac/Matrix.h>
#include <L0/Dirac/Gamma.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>
#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>
#include <L2/Input/FileReader.h>

/*
// for testing it can be useful to store smeared forward propagators
// #define __SAVE_SMEARED_PROPS__
*/

// comment this in if propagators have uniform temporal boundary contitions
// (e.g. the HMC inverter does this)
// for forward and sequential propagators separately
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__

// comment this if you don't want the threepoint function to be calculated
#define __CALCULATE_THREEPOINT__

// comment this if you don't want the twopoint function to be calculated
#define __CALCULATE_TWOPOINT__



int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);
  size_t L_tmp = 0;
  size_t T_tmp = 0;



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

  Input::FileReader reader("./simon_test_input.xml");

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
  int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
  size_t const timeslice_source = source_position[Base::idx_T] % T;
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
  size_t const timeslice_sink = (timeslice_source +  sourceSinkSeparation) % T;
  if (weave.isRoot())
    std::cout << "timeslice (sink) = " << timeslice_sink << std::endl;

  // make sure the boundary is not crossed by source-sink correlaton function
  size_t const timeslice_boundary = (timeslice_source + (T/2)) % T;
//   size_t const timeslice_boundary = (timeslice_sink + 1) % T;
  if (weave.isRoot())
    std::cout << "timeslice (boundary) = " << timeslice_boundary << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
  std::vector< std::string > const &seqPropFilesD(files[3]);
  std::vector< std::string > const &seqPropFilesU(files[4]);
/*
#ifdef __SAVE_SMEARED_PROPS__
  assert(files.size() >= 7);
  std::vector< std::string > const &SmearedPropFiles_d(files[5]);
  std::vector< std::string > const &SmearedPropFiles_u(files[6]);
#endif
*/
  
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


  Dirac::Gamma< 24 > g2g0;
  Dirac::Gamma< 245 > g2g0g5;
  Dirac::Gamma< 5 > g5;
  Dirac::Gamma< -1 > identity;
  std::complex< double > const factorP(+1.0, 0.0);

#ifdef __CALCULATE_TWOPOINT__

  if (weave.isRoot())
    std::cout << "\n calculating 2-point function \n" << std::endl;


  forwardProp_u.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  forwardProp_d.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  if (weave.isRoot())
    std::cout << "propagators and gauge field smeared successfully\n" << std::endl;

/*
#ifdef __SAVE_SMEARED_PROPS__
 Tool::IO::save(&forwardProp_u, SmearedPropFiles_u, Tool::IO::fileSCIDAC);
 Tool::IO::save(&forwardProp_d, SmearedPropFiles_d, Tool::IO::fileSCIDAC);
 if (weave.isRoot())
    std::cout << "smeared propagators saved successfully\n" << std::endl;
#endif
*/

  forwardProp_u.rotateToPhysicalBasis(true);
  forwardProp_d.rotateToPhysicalBasis(false);
  
  
  Core::BaryonCorrelator C2_P_11 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C2_N_11 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C2_P_22 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C2_N_22 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C2_P_12 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C2_N_12 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C2_P_21 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  
  Core::BaryonCorrelator C2_N_21 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  

  C2_P_11.prepareMomentumProjection(sourcePos);

  
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C2_P_11.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C2_N_11.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_2point_proton.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_2point_proton_projected.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron_projected.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C2_P_22.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C2_N_22.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_2point_proton.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_2point_proton_projected.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron_projected.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C2_P_21.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C2_N_21.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_2point_proton.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_2point_proton_projected.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron_projected.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C2_P_12.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C2_N_12.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_2point_proton.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_2point_proton_projected.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron_projected.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }

  weave.barrier();


#endif
// (if __CALCULATE_TWOPOINT__)

#ifdef __CALCULATE_THREEPOINT__

  Core::Propagator sequentialProp_u(L, T);
  Core::Propagator sequentialProp_d(L, T);

  Tool::IO::load(&sequentialProp_u, seqPropFilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark sequential propagator successfully loaded\n" << std::endl;

  Tool::IO::load(&sequentialProp_d, seqPropFilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark sequential propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__
  sequentialProp_u.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
  sequentialProp_d.changeBoundaryConditions_uniformToFixed(timeslice_sink, timeslice_boundary);
#endif


  sequentialProp_u.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  sequentialProp_d.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  if (weave.isRoot())
    std::cout << "sequential propagators smeared successfully\n" << std::endl;

  sequentialProp_u.rotateToPhysicalBasis(true);
  sequentialProp_d.rotateToPhysicalBasis(false);
  

  if (weave.isRoot())
    std::cout << "\n calculating 3-point function(s) \n" << std::endl;


  {  
  Core::BaryonCorrelator C3_P_DD_11 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, sequentialProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_UU_11 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, sequentialProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_DD_22 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, sequentialProp_d,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_UU_22 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, sequentialProp_u,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_DD_12 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, sequentialProp_d,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_UU_12 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, sequentialProp_u,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_DD_21 = Contract::nucleon_twopoint(forwardProp_u, forwardProp_u, sequentialProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_UU_21 = Contract::nucleon_twopoint(forwardProp_d, forwardProp_d, sequentialProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  

  C3_P_DD_11.prepareMomentumProjection(sourcePos);

  
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_DD_11.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_UU_11.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_dd.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_dd_projected.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu_projected.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_DD_22.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_UU_22.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_dd.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_dd_projected.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu_projected.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_DD_21.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_UU_21.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_dd.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_dd_projected.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu_projected.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_DD_12.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_UU_12.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_dd.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_dd_projected.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_uu_projected.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  
  } // dd current done

  {
  Core::BaryonCorrelator C3_P_UU_11 = Contract::nucleon_twopoint(sequentialProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_11 = Contract::nucleon_twopoint(sequentialProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_UU_22 = Contract::nucleon_twopoint(sequentialProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_22 = Contract::nucleon_twopoint(sequentialProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_UU_12 = Contract::nucleon_twopoint(sequentialProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_12 = Contract::nucleon_twopoint(sequentialProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_UU_21 = Contract::nucleon_twopoint(sequentialProp_u, forwardProp_u, forwardProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_21 = Contract::nucleon_twopoint(sequentialProp_d, forwardProp_d, forwardProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  

  Core::BaryonCorrelator C3_P_UU_11_tmp = Contract::nucleon_twopoint(forwardProp_u, sequentialProp_u, forwardProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_11_tmp = Contract::nucleon_twopoint(forwardProp_d, sequentialProp_d, forwardProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_UU_22_tmp = Contract::nucleon_twopoint(forwardProp_u, sequentialProp_u, forwardProp_d,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_22_tmp = Contract::nucleon_twopoint(forwardProp_d, sequentialProp_d, forwardProp_u,
                                                 factorP, 
                                                 g5   /* GammaSourceA */, 
                                                 g2g0 /* GammaSourceB */, 
                                                 g5   /* GammaSinkA */, 
                                                 g2g0 /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_UU_12_tmp = Contract::nucleon_twopoint(forwardProp_u, sequentialProp_u, forwardProp_d,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_12_tmp = Contract::nucleon_twopoint(forwardProp_d, sequentialProp_d, forwardProp_u,
                                                 factorP, 
                                                 g5       /* GammaSourceA */, 
                                                 g2g0     /* GammaSourceB */, 
                                                 identity /* GammaSinkA */, 
                                                 g2g0g5   /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_P_UU_21_tmp = Contract::nucleon_twopoint(forwardProp_u, sequentialProp_u, forwardProp_d,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  
  Core::BaryonCorrelator C3_N_DD_21_tmp = Contract::nucleon_twopoint(forwardProp_d, sequentialProp_d, forwardProp_u,
                                                 factorP, 
                                                 identity /* GammaSourceA */, 
                                                 g2g0g5   /* GammaSourceB */, 
                                                 g5       /* GammaSinkA */, 
                                                 g2g0     /* GammaSinkB */);
  
  C3_P_UU_11.prepareMomentumProjection(sourcePos);
  
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_UU_11.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_DD_11.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_tmp(C3_P_UU_11_tmp.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN_tmp(C3_N_DD_11_tmp.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_uu.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] += all_corrsP_tmp[I];
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] += all_corrsN_tmp[I];
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_uu_projected.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd_projected.11.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_UU_22.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_DD_22.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_tmp(C3_P_UU_22_tmp.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN_tmp(C3_N_DD_22_tmp.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_uu.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] += all_corrsP_tmp[I];
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] += all_corrsN_tmp[I];
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_uu_projected.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd_projected.22.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_UU_21.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_DD_21.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_tmp(C3_P_UU_21_tmp.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN_tmp(C3_N_DD_21_tmp.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_uu.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] += all_corrsP_tmp[I];
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] += all_corrsN_tmp[I];
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_uu_projected.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd_projected.21.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  {
    std::vector< Core::BaryonCorrelator > all_corrsP(C3_P_UU_12.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN(C3_N_DD_12.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_tmp(C3_P_UU_12_tmp.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsN_tmp(C3_N_DD_12_tmp.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_3point_proton_uu.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] += all_corrsP_tmp[I];
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] += all_corrsN_tmp[I];
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout =  new std::ofstream("output_3point_proton_uu_projected.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_3point_neutron_dd_projected.12.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    weave.barrier();
    delete fout;
  }
  
  } // uu current done

  
#endif
// (if __CALCULATE_THREEPOINT__)

  for(size_t I=0; I<momenta.size(); I++)
    delete [] momenta[I];
  momenta.clear();

  if (weave.isRoot())
    std::cout << "contractions performed and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
