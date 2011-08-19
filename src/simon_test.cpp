
#include <cstring>
#include <ctime>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Ahmidas.h>
#include <L0/Base/Weave.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Matrix.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>
#include <L2/Input/FileReader.h>


int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  time_t start,end;
  double dif;

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("simon_test_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::vector< int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);

  size_t const L(L_tmp);
  size_t const T(T_tmp);
  Base::Weave weave(L, T);

  std::ofstream *fout = NULL;

  if (weave.isRoot())
  std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu    = floats["mu"];
  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = floats["APE_steps"];

  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = floats["Jac_steps"];

  size_t const timeslice_source = (positions[0])[Base::idx_T];
  size_t const timeslice_boundary((timeslice_source + T/2) % T);
  size_t const timeslice_stochSource = (timeslice_source + size_t(floats["sourceSinkSeparation"])) % T;
  assert(timeslice_stochSource != timeslice_source);
  assert(timeslice_source < T);
  assert(timeslice_boundary < T);
  assert(timeslice_stochSource < T);
  size_t const timeslice_sink = timeslice_stochSource;

  if (weave.isRoot())
  {
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;
    std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
    std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
    std::cout << "timeslice (stochastic wall source) = " << timeslice_stochSource << std::endl;
  }


  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);
  std::vector< std::string > const &stochasticPropFilesD(files[3]);
  std::vector< std::string > const &stochasticPropFilesU(files[4]);
  std::vector< std::string > const &stochasticSourceFilesD(files[5]);
  std::vector< std::string > const &stochasticSourceFilesU(files[6]);


  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Core::Field< QCD::Gauge > gauge_field(L, T);
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;


  Core::Propagator uProp = Core::Propagator(L, T);

  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC);

  if (weave.isRoot())
    std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp = Core::Propagator(L, T);

  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC);

  if (weave.isRoot())
    std::cout << "d quark propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 12 > stochastic_dProp(L, T);
  Core::StochasticPropagator< 12 > stochastic_uProp(L, T);
  Core::StochasticSource< 12 > stochasticSourceD(L, T);
  Core::StochasticSource< 12 > stochasticSourceU(L, T);


  // the load function with last parameter (size_t) precision is sort of a quick and dirty version
  // but it also reads files without a proper header
  Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochastic_dProp),
                 stochasticPropFilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "stochastic d quark propagator successfully loaded\n" << std::endl;

  Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochastic_uProp),
                 stochasticPropFilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "stochastic u quark propagator successfully loaded\n" << std::endl;

  Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochasticSourceD),
                 stochasticSourceFilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "stochastic source (d) successfully loaded\n" << std::endl;

  Tool::IO::load(dynamic_cast< Core::Propagator *> (&stochasticSourceU),
                 stochasticSourceFilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "stochastic source (u) successfully loaded\n" << std::endl;

  // change boundary conditions

  dProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  uProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);

  Core::Propagator uProp_smeared(uProp);
  Core::Propagator dProp_smeared(dProp);

  
  Smear::APE APE_tool(APE_alpha);
  Smear::Jacobi Jacobi_tool(Jac_alpha);

  // we still need the unsmeared gauge field for the contractions, therefore we make a copy here
  Core::Field< QCD::Gauge > gauge_field_APE(gauge_field);
  APE_tool.smear(gauge_field_APE, APE_iterations);

  uProp_smeared.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_APE);
  dProp_smeared.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_APE);
  

  if (weave.isRoot())
    std::cout << "forward propagators smeared successfully\n" << std::endl;
  

  stochastic_dProp.changeBoundaryConditions_uniformToFixed(timeslice_stochSource, timeslice_boundary);
  stochastic_uProp.changeBoundaryConditions_uniformToFixed(timeslice_stochSource, timeslice_boundary);

  std::vector< Base::Operator > my_operators;
  my_operators.push_back(Base::op_GAMMA_15);
  my_operators.push_back(Base::op_GAMMA_25);
  my_operators.push_back(Base::op_GAMMA_35);
//   my_operators.push_back(Base::op_O44);
//   my_operators.push_back(Base::op_O11);
//   my_operators.push_back(Base::op_O22);
//   my_operators.push_back(Base::op_O33);
  // my_operators.push_back(Base::op_UNITY);

  // ****************************************************************************************************
  // ****************************************************************************************************
  // ****************************************************************************************************
  
  std::vector< Core::BaryonCorrelator > p3p;
  
  
  // STOCHASTIC TWOPOINT VERSION 1

/*
  if (weave.isRoot())
    std::cout << "starting contractions version 1" << std::endl;
  
  weave.barrier();
  time (&start);
                

  p3p = Contract::proton_threepoint_stochastic_naive(uProp_smeared, dProp_smeared,
                                                          uProp, dProp,
                                                          stochastic_uProp, stochastic_dProp,
                                                          stochasticSourceU, stochasticSourceD,
                                                          NULL, //&gauge_field,
                                                          my_operators,
                                                          timeslice_source, timeslice_sink);
  
  weave.barrier();
  time (&end);
  dif = difftime (end,start);
  if (weave.isRoot())
    std::cout << "done in " << dif << "seconds\n" << std::endl;

  p3p[0] *= Base::proj_1_PLUS_TM;
  p3p[1] *= Base::proj_1_PLUS_TM;
  p3p[2] *= Base::proj_2_PLUS_TM;
  p3p[3] *= Base::proj_2_PLUS_TM;
  p3p[4] *= Base::proj_3_PLUS_TM;
  p3p[5] *= Base::proj_3_PLUS_TM;


  if (weave.isRoot())
  {
    fout = new std::ofstream("output_3point_proton_stochastic_naive_uu.dat");
    for (size_t i=0; i<p3p.size(); i+=2)
    {
      p3p[i].setOffset(timeslice_source);
      if (weave.isRoot())
      {
        *fout << p3p[i] << std::endl;
      }
    }
    fout->close();
    fout->open("output_3point_proton_stochastic_naive_dd.dat");
    for (size_t i=1; i<p3p.size(); i+=2)
    {
      p3p[i].setOffset(timeslice_source);
      if (weave.isRoot())
      {
        *fout << p3p[i] << std::endl;
      }
    }
    fout->close();
  }
  weave.barrier();
  p3p.clear();
*/

  // ****************************************************************************************************
  // ****************************************************************************************************
  // ****************************************************************************************************

  // STOCHASTIC TWOPOINT VERSION 2

  if (weave.isRoot())
    std::cout << "starting contractions version 2\n" << std::endl;
  
  {
  weave.barrier();
  time (&start);


  my_operators.push_back(Base::op_O44);
  my_operators.push_back(Base::op_O11);
  my_operators.push_back(Base::op_O22);
  my_operators.push_back(Base::op_O33);
  std::vector< Base::BaryonPropagatorProjector > my_projectors_u;
  my_projectors_u.push_back(Base::proj_1_PLUS_TM);
  my_projectors_u.push_back(Base::proj_2_PLUS_TM);
  my_projectors_u.push_back(Base::proj_3_PLUS_TM);
  my_projectors_u.push_back(Base::proj_PARITY_PLUS_TM);
  my_projectors_u.push_back(Base::proj_PARITY_PLUS_TM);
  my_projectors_u.push_back(Base::proj_PARITY_PLUS_TM);
  my_projectors_u.push_back(Base::proj_PARITY_PLUS_TM);
  std::vector< Base::BaryonPropagatorProjector > my_projectors_d(my_projectors_u);
//   my_projectors_d.push_back(Base::proj_1_MINUS_TM);
//   my_projectors_d.push_back(Base::proj_2_MINUS_TM);
//   my_projectors_d.push_back(Base::proj_3_MINUS_TM);

  p3p = Contract::proton_threepoint_stochastic_non_local(
//   p3p = Contract::proton_threepoint_stochastic(
                                               uProp_smeared, dProp_smeared, uProp, dProp,
                                               stochastic_uProp, stochastic_dProp,
                                               stochasticSourceU, stochasticSourceD,
                                               &gauge_field,
                                               timeslice_source, timeslice_stochSource,
                                               my_operators, my_projectors_u, my_projectors_d);

  // FIX THIS: problem with Gamma multiplication
  p3p[2*1]   *= -1.0;
  p3p[2*1+1] *= -1.0;
  p3p[2*5]   *= -1.0;
  p3p[2*5+1] *= -1.0;



  weave.barrier();
  time (&end);
  dif = difftime (end,start);
  if (weave.isRoot())
    std::cout << "done in " << dif << "seconds\n" << std::endl;

  if (weave.isRoot())
  {
    fout = new std::ofstream("output_3point_proton_stochastic_uu.dat");
    for (size_t i=0; i<p3p.size(); i+=2)
    {
      p3p[i].setOffset(timeslice_source);
      if (weave.isRoot())
      {
        *fout << p3p[i] << std::endl;
      }
    }
    fout->close();
    fout->open("output_3point_proton_stochastic_dd.dat");
    for (size_t i=1; i<p3p.size(); i+=2)
    {
      p3p[i].setOffset(timeslice_source);
      if (weave.isRoot())
      {
        *fout << p3p[i] << std::endl;
      }
    }
    fout->close();
  }
  p3p.clear();
  }

/*
  if (weave.isRoot())
    std::cout << "starting contractions version 2, non-local\n" << std::endl;
{
  weave.barrier();
  time (&start);

//   assert (my_operators.size() >= 5);
  std::vector< Base::BaryonPropagatorProjector > my_projectors_u;
  my_projectors_u.push_back(Base::proj_1_PLUS_TM);
  my_projectors_u.push_back(Base::proj_2_PLUS_TM);
  my_projectors_u.push_back(Base::proj_3_PLUS_TM);
  std::vector< Base::BaryonPropagatorProjector > my_projectors_d(my_projectors_u);
//   my_projectors_d.push_back(Base::proj_1_MINUS_TM);
//   my_projectors_d.push_back(Base::proj_2_MINUS_TM);
//   my_projectors_d.push_back(Base::proj_3_MINUS_TM);

  p3p = Contract::proton_threepoint_stochastic_non-local(uProp_smeared, dProp_smeared, uProp, dProp,
                                               stochastic_uProp, stochastic_dProp,
                                               stochasticSourceU, stochasticSourceD,
                                               timeslice_source, timeslice_stochSource,
                                               my_operators, my_projectors_u, my_projectors_d);

  weave.barrier();
  time (&end);
  dif = difftime (end,start);
  if (weave.isRoot())
    std::cout << "done in " << dif << "seconds\n" << std::endl;

  if (weave.isRoot())
  {
    fout = new std::ofstream("output_3point_proton_stochastic_uu.dat");
    for (size_t i=0; i<p3p.size(); i+=2)
    {
      p3p[i].setOffset(timeslice_source);
      if (weave.isRoot())
      {
        *fout << p3p[i] << std::endl;
      }
    }
    fout->close();
    fout->open("output_3point_proton_stochastic_dd.dat");
    for (size_t i=1; i<p3p.size(); i+=2)
    {
      p3p[i].setOffset(timeslice_source);
      if (weave.isRoot())
      {
        *fout << p3p[i] << std::endl;
      }
    }
    fout->close();
  }
  p3p.clear();
  }
*/

  // ****************************************************************************************************
  // ****************************************************************************************************
  // ****************************************************************************************************

  // 2-POINT FUNCTION

  if (weave.isRoot())
    std::cout << "starting 2-point contractions\n" << std::endl;

 {

  // now that the threepoint contractions are done, we can smear the stochastic propagators
  stochastic_uProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_APE);
  stochastic_dProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_APE);

  size_t const * const source_position = positions[0];

  // stochastic estimates of the real propagator
  // note flavour change here because we 'revert' the twisted mass propagator
  Core::Propagator dProp_stochEst((stochasticSourceU.createStochasticPropagator_fixedSink(stochastic_uProp, source_position)).revert());
  Core::Propagator uProp_stochEst((stochasticSourceD.createStochasticPropagator_fixedSink(stochastic_dProp, source_position)).revert());


  uProp_stochEst.rotateToPhysicalBasis(true);
  dProp_stochEst.rotateToPhysicalBasis(false);

  // int const * const dummy_momentum({0, 0, 0});

  uProp_smeared.rotateToPhysicalBasis(true);
  dProp_smeared.rotateToPhysicalBasis(false);

  // STOCHASTIC and EXACT 2-POINT FUNCTION

  Core::BaryonCorrelator C2_P = Contract::proton_twopoint(uProp_smeared, uProp_smeared, dProp_smeared, Base::proj_NO_PROJECTOR);

  // d line stochastically estimated
  Core::BaryonCorrelator C2_P_stochD     = Contract::proton_twopoint(uProp_smeared, uProp_smeared, dProp_stochEst, Base::proj_NO_PROJECTOR);
  C2_P_stochD.deleteField();

  // one u line stochastically estimated (symmetrized)
  Core::BaryonCorrelator C2_P_stochU     = Contract::proton_twopoint(uProp_stochEst, uProp_smeared, dProp_smeared, Base::proj_NO_PROJECTOR);
  Core::BaryonCorrelator C2_P_stochU_tmp = Contract::proton_twopoint(uProp_smeared, uProp_stochEst, dProp_smeared, Base::proj_NO_PROJECTOR);
  C2_P_stochU_tmp.deleteField();
  C2_P_stochU.deleteField();
  C2_P_stochU += C2_P_stochU_tmp;
  C2_P_stochU *= 0.5;


  C2_P_stochD.setOffset(timeslice_source);
  C2_P_stochD *= Base::proj_PARITY_PLUS_STD;
  if (weave.isRoot())
  {
    fout = new std::ofstream("output_2point_proton_stochD_projected.dat");
    *fout << C2_P_stochD << std::endl;
    fout->close();
  }

  C2_P_stochU.setOffset(timeslice_source);
  C2_P_stochU *= Base::proj_PARITY_PLUS_STD;
  if (weave.isRoot())
  {
    fout = new std::ofstream("output_2point_proton_stochU_projected.dat");
    *fout << C2_P_stochU << std::endl;
    fout->close();
  }

  C2_P.setOffset(timeslice_source);
  C2_P *= Base::proj_PARITY_PLUS_STD;
  if (weave.isRoot())
  {
    fout = new std::ofstream("output_2point_proton_projected.dat");
    *fout << C2_P << std::endl;
    fout->close();
  }

/*
  Core::BaryonCorrelator C2_N = Contract::proton_twopoint(dProp, dProp, uProp, Base::proj_NO_PROJECTOR);

  int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
  C2_P.prepareMomentumProjection(sourcePos);


  std::vector< Core::BaryonCorrelator > all_corrsP(C2_P.momentumProjection(momenta));

  std::vector< Core::BaryonCorrelator > all_corrsN(C2_N.momentumProjection(momenta));

  std::ofstream *fout = NULL;
  if (weave.isRoot())
  {
    fout = new std::ofstream("output_2point_proton.dat");
    for(size_t I=0; I<momenta.size(); I++)
    {
      all_corrsP[I].setOffset(timeslice_source);
      all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
    }
    fout->close();
    fout->open("output_2point_neutron.dat");
    for(size_t I=0; I<momenta.size(); I++)
    {
      all_corrsN[I].setOffset(timeslice_source);
      all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
    }
    fout->close();
    fout->open("output_2point_proton_projected.dat");
    for(size_t I=0; I<momenta.size(); I++)
    {
      all_corrsP[I].setOffset(timeslice_source);
      all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
      all_corrsP[I].printWithMomentum(*fout, momenta[I]);
    }
    fout->close();
    fout->open("output_2point_neutron_projected.dat");
    for(size_t I=0; I<momenta.size(); I++)
    {
      all_corrsN[I].setOffset(timeslice_source);
      all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
      all_corrsN[I].printWithMomentum(*fout, momenta[I]);
    }
    fout->close();
  }
*/
  }

  weave.barrier();
  if (weave.isRoot())
    std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  return EXIT_SUCCESS;
}
