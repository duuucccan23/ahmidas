
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

  if (weave.isRoot())
  std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu    = floats["mu"];
  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = floats["APE_steps"];

  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = floats["Jac_steps"];

  size_t const * const source_position = positions[0];
  size_t const timeslice_source = (positions[0])[Base::idx_T];
  size_t const timeslice_boundary((timeslice_source + T/2) % T);
  size_t const timeslice_stochSource = (timeslice_source + size_t(floats["sourceSinkSeparation"])) % T;
  assert(timeslice_stochSource != timeslice_source);
  assert(timeslice_source < T);
  assert(timeslice_boundary < T);
  assert(timeslice_stochSource < T);
  //size_t const timeslice_sink = timeslice_stochSource;

  if (weave.isRoot())
  {
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;
    std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
    std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;
    std::cout << "timeslice (stochastic wall source) = " << timeslice_stochSource << std::endl;
  }

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

  // ****************************************************************************************************
  // ****************************************************************************************************
  // ****************************************************************************************************

  // ****************************************************************************************************
  // ****************************************************************************************************
  // ****************************************************************************************************

  // STOCHASTIC THREEPOINT


  // ********************************************************************************************
  // local currents *****************************************************************************
  // ********************************************************************************************

  if (weave.isRoot())
    std::cout << "\n local currents (proton) ... \n" << std::endl;

  {
    weave.barrier();
    time (&start);

    std::vector< Core::BaryonCorrelator > C3p_dummy;

    std::vector< Core::BaryonCorrelator > C3p_proj0_local
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_PARITY_PLUS_TM, C3p_dummy);
    std::vector< Core::BaryonCorrelator > C3p_proj1_local
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_1_PLUS_TM, C3p_dummy);
    std::vector< Core::BaryonCorrelator > C3p_proj2_local
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_2_PLUS_TM, C3p_dummy);
    std::vector< Core::BaryonCorrelator > C3p_proj3_local
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_3_PLUS_TM, C3p_dummy);

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C3p_proj0_local[0].prepareMomentumProjection(sourcePos);

    std::ofstream *fout_proj0_uu = NULL;
    std::ofstream *fout_proj0_dd = NULL;
    std::ofstream *fout_proj1_uu = NULL;
    std::ofstream *fout_proj1_dd = NULL;
    std::ofstream *fout_proj2_uu = NULL;
    std::ofstream *fout_proj2_dd = NULL;
    std::ofstream *fout_proj3_uu = NULL;
    std::ofstream *fout_proj3_dd = NULL;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_proton_stochastic_local_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_proton_stochastic_local_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_proton_stochastic_local_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_proton_stochastic_local_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_proton_stochastic_local_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_proton_stochastic_local_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_proton_stochastic_local_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_proton_stochastic_local_proj3_dd.dat");
    }

    assert(C3p_proj0_local.size() == 32);

    for(size_t i=0; i<C3p_proj0_local.size(); i++)
    {
      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_local[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_local[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_local[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_local[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }


  if (weave.isRoot())
    std::cout << "\n local currents (neutron) ... \n" << std::endl;

  {
    weave.barrier();
    time (&start);

    std::vector< Core::BaryonCorrelator > C3n_dummy;

    std::vector< Core::BaryonCorrelator > C3n_proj0_local
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_PARITY_PLUS_TM, C3n_dummy);
    std::vector< Core::BaryonCorrelator > C3n_proj1_local
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_1_PLUS_TM, C3n_dummy);
    std::vector< Core::BaryonCorrelator > C3n_proj2_local
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_2_PLUS_TM, C3n_dummy);
    std::vector< Core::BaryonCorrelator > C3n_proj3_local
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "local", Base::proj_3_PLUS_TM, C3n_dummy);

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C3n_proj0_local[0].prepareMomentumProjection(sourcePos);

    std::ofstream *fout_proj0_uu = NULL;
    std::ofstream *fout_proj0_dd = NULL;
    std::ofstream *fout_proj1_uu = NULL;
    std::ofstream *fout_proj1_dd = NULL;
    std::ofstream *fout_proj2_uu = NULL;
    std::ofstream *fout_proj2_dd = NULL;
    std::ofstream *fout_proj3_uu = NULL;
    std::ofstream *fout_proj3_dd = NULL;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_neutron_stochastic_local_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_neutron_stochastic_local_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_neutron_stochastic_local_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_neutron_stochastic_local_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_neutron_stochastic_local_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_neutron_stochastic_local_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_neutron_stochastic_local_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_neutron_stochastic_local_proj3_dd.dat");
    }

    assert(C3n_proj0_local.size() == 32);

    for(size_t i=0; i<C3n_proj0_local.size(); i++)
    {

      // ***************************************************************
      // *** nota bene: u and d currents are swapped for the neutron ***
      // *** therefore we have to correct for tau_3 factors ************
      // ***************************************************************
      switch (i)
      {
        case  0:
        case  1:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
        case 17:
        case 26:
        case 27:
        case 28:
        case 29:
        case 30:
        case 31:
          C3n_proj0_local[i] *= -1.0;
          C3n_proj1_local[i] *= -1.0;
          C3n_proj2_local[i] *= -1.0;
          C3n_proj3_local[i] *= -1.0;
          break;
        default:
          // do nothing
          break;
      }

      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3n_proj0_local[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3n_proj1_local[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3n_proj2_local[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3n_proj3_local[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          // since u and d currents are swapped for the neutron the d current comes first
          // which is contrary to the proton threepoint
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }


  // ********************************************************************************************
  // conserved currents *************************************************************************
  // ********************************************************************************************

  if (weave.isRoot())
    std::cout << "\n conserved currents (proton) ... \n" << std::endl;

  {
    weave.barrier();
    time (&start);

    std::vector< Core::BaryonCorrelator > C3p_proj0_conservedA;
    std::vector< Core::BaryonCorrelator > C3p_proj1_conservedA;
    std::vector< Core::BaryonCorrelator > C3p_proj2_conservedA;
    std::vector< Core::BaryonCorrelator > C3p_proj3_conservedA;

    std::vector< Core::BaryonCorrelator > C3p_proj0_conserved
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_PARITY_PLUS_TM, C3p_proj0_conservedA);
    std::vector< Core::BaryonCorrelator > C3p_proj1_conserved
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_1_PLUS_TM, C3p_proj1_conservedA);
    std::vector< Core::BaryonCorrelator > C3p_proj2_conserved
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_2_PLUS_TM, C3p_proj2_conservedA);
    std::vector< Core::BaryonCorrelator > C3p_proj3_conserved
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_3_PLUS_TM, C3p_proj3_conservedA);

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C3p_proj0_conserved[0].prepareMomentumProjection(sourcePos);

    std::ofstream *fout_proj0_uu = NULL;
    std::ofstream *fout_proj0_dd = NULL;
    std::ofstream *fout_proj1_uu = NULL;
    std::ofstream *fout_proj1_dd = NULL;
    std::ofstream *fout_proj2_uu = NULL;
    std::ofstream *fout_proj2_dd = NULL;
    std::ofstream *fout_proj3_uu = NULL;
    std::ofstream *fout_proj3_dd = NULL;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_vector_proj3_dd.dat");
    }

    assert(C3p_proj0_conserved.size() == 8);

    for(size_t i=0; i<C3p_proj0_conserved.size(); i++)
    {
      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_conserved[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_conserved[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_conserved[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_conserved[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_proton_stochastic_conserved_axial_proj3_dd.dat");
    }

    assert(C3p_proj0_conservedA.size() == 8);

    for(size_t i=0; i<C3p_proj0_conserved.size(); i++)
    {
      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_conservedA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_conservedA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_conservedA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_conservedA[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }


  if (weave.isRoot())
    std::cout << "\n conserved currents (neutron) ... \n" << std::endl;

  {
    weave.barrier();
    time (&start);

    std::vector< Core::BaryonCorrelator > C3n_proj0_conservedA;
    std::vector< Core::BaryonCorrelator > C3n_proj1_conservedA;
    std::vector< Core::BaryonCorrelator > C3n_proj2_conservedA;
    std::vector< Core::BaryonCorrelator > C3n_proj3_conservedA;

    std::vector< Core::BaryonCorrelator > C3n_proj0_conserved
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_PARITY_PLUS_TM, C3n_proj0_conservedA);
    std::vector< Core::BaryonCorrelator > C3n_proj1_conserved
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_1_PLUS_TM, C3n_proj1_conservedA);
    std::vector< Core::BaryonCorrelator > C3n_proj2_conserved
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_2_PLUS_TM, C3n_proj2_conservedA);
    std::vector< Core::BaryonCorrelator > C3n_proj3_conserved
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "noether", Base::proj_3_PLUS_TM, C3n_proj3_conservedA);

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C3n_proj0_conserved[0].prepareMomentumProjection(sourcePos);

    std::ofstream *fout_proj0_uu = NULL;
    std::ofstream *fout_proj0_dd = NULL;
    std::ofstream *fout_proj1_uu = NULL;
    std::ofstream *fout_proj1_dd = NULL;
    std::ofstream *fout_proj2_uu = NULL;
    std::ofstream *fout_proj2_dd = NULL;
    std::ofstream *fout_proj3_uu = NULL;
    std::ofstream *fout_proj3_dd = NULL;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_vector_proj3_dd.dat");
    }

    assert(C3n_proj0_conserved.size() == 8);

    for(size_t i=0; i<C3n_proj0_conserved.size(); i++)
    {

      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3n_proj0_conserved[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3n_proj1_conserved[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3n_proj2_conserved[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3n_proj3_conserved[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          // since u and d currents are swapped for the neutron the d current comes first
          // which is contrary to the proton threepoint
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_neutron_stochastic_conserved_axial_proj3_dd.dat");
    }

    assert(C3n_proj0_conservedA.size() == 8);

    for(size_t i=0; i<C3n_proj0_conserved.size(); i++)
    {

      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3n_proj0_conservedA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3n_proj1_conservedA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3n_proj2_conservedA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3n_proj3_conservedA[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          // since u and d currents are swapped for the neutron the d current comes first
          // which is contrary to the proton threepoint
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }


  // ********************************************************************************************
  // derivative operators ***********************************************************************
  // ********************************************************************************************

  if (weave.isRoot())
    std::cout << "\n 1D currents (proton) ... \n" << std::endl;

  {
    weave.barrier();
    time (&start);

    std::vector< Core::BaryonCorrelator > C3p_proj0_1DA;
    std::vector< Core::BaryonCorrelator > C3p_proj1_1DA;
    std::vector< Core::BaryonCorrelator > C3p_proj2_1DA;
    std::vector< Core::BaryonCorrelator > C3p_proj3_1DA;

    std::vector< Core::BaryonCorrelator > C3p_proj0_1D
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_PARITY_PLUS_TM, C3p_proj0_1DA);
    std::vector< Core::BaryonCorrelator > C3p_proj1_1D
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_1_PLUS_TM, C3p_proj1_1DA);
    std::vector< Core::BaryonCorrelator > C3p_proj2_1D
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_2_PLUS_TM, C3p_proj2_1DA);
    std::vector< Core::BaryonCorrelator > C3p_proj3_1D
       = Contract::proton_threepoint_stochastic(uProp_smeared, dProp_smeared, uProp, dProp,
                                                stochastic_uProp, stochastic_dProp,
                                                stochasticSourceU, stochasticSourceD,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_3_PLUS_TM, C3p_proj3_1DA);

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C3p_proj0_1D[0].prepareMomentumProjection(sourcePos);

    std::ofstream *fout_proj0_uu = NULL;
    std::ofstream *fout_proj0_dd = NULL;
    std::ofstream *fout_proj1_uu = NULL;
    std::ofstream *fout_proj1_dd = NULL;
    std::ofstream *fout_proj2_uu = NULL;
    std::ofstream *fout_proj2_dd = NULL;
    std::ofstream *fout_proj3_uu = NULL;
    std::ofstream *fout_proj3_dd = NULL;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_proton_stochastic_1D_unpolarized_proj3_dd.dat");
    }

    assert(C3p_proj0_1D.size() == 32);

    for(size_t i=0; i<C3p_proj0_1D.size(); i++)
    {
      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_1D[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_1D[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_1D[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_1D[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_proton_stochastic_1D_polarized_proj3_dd.dat");
    }

    assert(C3p_proj0_1DA.size() == 32);

    for(size_t i=0; i<C3p_proj0_1D.size(); i++)
    {
      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3p_proj0_1DA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3p_proj1_1DA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3p_proj2_1DA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3p_proj3_1DA[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }


  if (weave.isRoot())
    std::cout << "\n 1D currents (neutron) ... \n" << std::endl;

  {
    weave.barrier();
    time (&start);

    std::vector< Core::BaryonCorrelator > C3n_proj0_1DA;
    std::vector< Core::BaryonCorrelator > C3n_proj1_1DA;
    std::vector< Core::BaryonCorrelator > C3n_proj2_1DA;
    std::vector< Core::BaryonCorrelator > C3n_proj3_1DA;

    std::vector< Core::BaryonCorrelator > C3n_proj0_1D
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_PARITY_PLUS_TM, C3n_proj0_1DA);
    std::vector< Core::BaryonCorrelator > C3n_proj1_1D
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_1_PLUS_TM, C3n_proj1_1DA);
    std::vector< Core::BaryonCorrelator > C3n_proj2_1D
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_2_PLUS_TM, C3n_proj2_1DA);
    std::vector< Core::BaryonCorrelator > C3n_proj3_1D
       = Contract::proton_threepoint_stochastic(dProp_smeared, uProp_smeared, dProp, uProp,
                                                stochastic_dProp, stochastic_uProp,
                                                stochasticSourceD, stochasticSourceU,
                                                &gauge_field,
                                                timeslice_source, timeslice_stochSource,
                                                "1D", Base::proj_3_PLUS_TM, C3n_proj3_1DA);

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C3n_proj0_1D[0].prepareMomentumProjection(sourcePos);

    std::ofstream *fout_proj0_uu = NULL;
    std::ofstream *fout_proj0_dd = NULL;
    std::ofstream *fout_proj1_uu = NULL;
    std::ofstream *fout_proj1_dd = NULL;
    std::ofstream *fout_proj2_uu = NULL;
    std::ofstream *fout_proj2_dd = NULL;
    std::ofstream *fout_proj3_uu = NULL;
    std::ofstream *fout_proj3_dd = NULL;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_unpolarized_proj3_dd.dat");
    }

    assert(C3n_proj0_1D.size() == 32);

    for(size_t i=0; i<C3n_proj0_1D.size(); i++)
    {

      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3n_proj0_1D[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3n_proj1_1D[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3n_proj2_1D[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3n_proj3_1D[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          // since u and d currents are swapped for the neutron the d current comes first
          // which is contrary to the proton threepoint
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    if (weave.isRoot())
    {
      fout_proj0_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj0_uu.dat");
      fout_proj0_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj0_dd.dat");
      fout_proj1_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj1_uu.dat");
      fout_proj1_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj1_dd.dat");
      fout_proj2_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj2_uu.dat");
      fout_proj2_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj2_dd.dat");
      fout_proj3_uu =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj3_uu.dat");
      fout_proj3_dd =  new std::ofstream("output_3point_neutron_stochastic_1D_polarized_proj3_dd.dat");
    }

    assert(C3n_proj0_1DA.size() == 32);

    for(size_t i=0; i<C3n_proj0_1D.size(); i++)
    {

      std::vector< Core::BaryonCorrelator > all_corrs_proj0(C3n_proj0_1DA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj1(C3n_proj1_1DA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj2(C3n_proj2_1DA[i].momentumProjection(momenta));
      std::vector< Core::BaryonCorrelator > all_corrs_proj3(C3n_proj3_1DA[i].momentumProjection(momenta));

      for(size_t I=0; I<momenta.size(); I++)
      {
        std::ostringstream oss;
        oss << std::setw(3) << i/2 << " " << std::flush;
        std::string const prefix(oss.str());
        all_corrs_proj0[I].setOffset(timeslice_source);
        all_corrs_proj1[I].setOffset(timeslice_source);
        all_corrs_proj2[I].setOffset(timeslice_source);
        all_corrs_proj3[I].setOffset(timeslice_source);
        if (weave.isRoot())
        {
          // since u and d currents are swapped for the neutron the d current comes first
          // which is contrary to the proton threepoint
          if(i%2 == 0)
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_dd, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_dd, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_dd, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_dd, momenta[I], prefix);
          }
          else
          {
            all_corrs_proj0[I].printWithMomentum(*fout_proj0_uu, momenta[I], prefix);
            all_corrs_proj1[I].printWithMomentum(*fout_proj1_uu, momenta[I], prefix);
            all_corrs_proj2[I].printWithMomentum(*fout_proj2_uu, momenta[I], prefix);
            all_corrs_proj3[I].printWithMomentum(*fout_proj3_uu, momenta[I], prefix);
          }
        }
      }
      weave.barrier();
    }

    if (weave.isRoot())
    {
      fout_proj0_uu->close();
      fout_proj0_dd->close();
      fout_proj1_uu->close();
      fout_proj1_dd->close();
      fout_proj2_uu->close();
      fout_proj2_dd->close();
      fout_proj3_uu->close();
      fout_proj3_dd->close();
    }

    delete fout_proj0_uu;
    delete fout_proj0_dd;
    delete fout_proj1_uu;
    delete fout_proj1_dd;
    delete fout_proj2_uu;
    delete fout_proj2_dd;
    delete fout_proj3_uu;
    delete fout_proj3_dd;

    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }


  // ****************************************************************************************************
  // ****************************************************************************************************
  // ****************************************************************************************************

  // 2-POINT FUNCTION


  // now that the threepoint contractions are done, we can smear the stochastic propagators
  stochastic_uProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_APE);
  stochastic_dProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_APE);

  uProp_smeared.rotateToPhysicalBasis(true);
  dProp_smeared.rotateToPhysicalBasis(false);

  // stochastic estimates of the real propagator
  // note flavour change here because we 'revert' the twisted mass propagator
  Core::Propagator dProp_stochEst((stochasticSourceU.createStochasticPropagator_fixedSink(stochastic_uProp, source_position)).revert());
  Core::Propagator uProp_stochEst((stochasticSourceD.createStochasticPropagator_fixedSink(stochastic_dProp, source_position)).revert());

  uProp_stochEst.rotateToPhysicalBasis(true);
  dProp_stochEst.rotateToPhysicalBasis(false);

  {
    if (weave.isRoot())
      std::cout << "starting 2-point contractions (proton)\n" << std::endl;

    weave.barrier();
    time (&start);

    // STOCHASTIC and EXACT 2-POINT FUNCTION

    Core::BaryonCorrelator C2_P = Contract::proton_twopoint(uProp_smeared, uProp_smeared, dProp_smeared, Base::proj_NO_PROJECTOR);

    // d line stochastically estimated
    Core::BaryonCorrelator C2_P_stochD     = Contract::proton_twopoint(uProp_smeared, uProp_smeared, dProp_stochEst, Base::proj_NO_PROJECTOR);

    // one u line stochastically estimated (symmetrized)
    Core::BaryonCorrelator C2_P_stochU     = Contract::proton_twopoint(uProp_stochEst, uProp_smeared, dProp_smeared, Base::proj_NO_PROJECTOR);
    Core::BaryonCorrelator C2_P_stochU_tmp = Contract::proton_twopoint(uProp_smeared, uProp_stochEst, dProp_smeared, Base::proj_NO_PROJECTOR);
    C2_P_stochU += C2_P_stochU_tmp;
    C2_P_stochU *= 0.5;

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C2_P.prepareMomentumProjection(sourcePos);

    std::vector< Core::BaryonCorrelator > all_corrsP(C2_P.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_stochD(C2_P_stochD.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_stochU(C2_P_stochU.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_2point_proton.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout->open("output_2point_proton_stochD.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP_stochD[I].setOffset(timeslice_source);
        all_corrsP_stochD[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout->open("output_2point_proton_stochU.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP_stochU[I].setOffset(timeslice_source);
        all_corrsP_stochU[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      delete fout;
    }

    C2_P_stochD.setOffset(timeslice_source);
    C2_P_stochD *= Base::proj_PARITY_PLUS_STD;
    if (weave.isRoot())
    {
      fout = new std::ofstream("output_2point_proton_stochD_projected.dat");
      *fout << C2_P_stochD << std::endl;
      fout->close();
      delete fout;
    }
    C2_P_stochU.setOffset(timeslice_source);
    C2_P_stochU *= Base::proj_PARITY_PLUS_STD;
    if (weave.isRoot())
    {
      fout = new std::ofstream("output_2point_proton_stochU_projected.dat");
      *fout << C2_P_stochU << std::endl;
      fout->close();
      delete fout;
    }
    C2_P.setOffset(timeslice_source);
    C2_P *= Base::proj_PARITY_PLUS_STD;
    if (weave.isRoot())
    {
      fout = new std::ofstream("output_2point_proton_projected.dat");
      *fout << C2_P << std::endl;
      fout->close();
      delete fout;
    }
    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }

  {
    if (weave.isRoot())
      std::cout << "starting 2-point contractions (neutron)\n" << std::endl;

    weave.barrier();
    time (&start);

    // STOCHASTIC and EXACT 2-POINT FUNCTION

    Core::BaryonCorrelator C2_N = Contract::proton_twopoint(dProp_smeared, dProp_smeared, uProp_smeared, Base::proj_NO_PROJECTOR);

    // d line stochastically estimated
    Core::BaryonCorrelator C2_N_stochU     = Contract::proton_twopoint(dProp_smeared, dProp_smeared, uProp_stochEst, Base::proj_NO_PROJECTOR);

    // one u line stochastically estimated (symmetrized)
    Core::BaryonCorrelator C2_N_stochD     = Contract::proton_twopoint(dProp_stochEst, dProp_smeared, uProp_smeared, Base::proj_NO_PROJECTOR);
    Core::BaryonCorrelator C2_N_stochD_tmp = Contract::proton_twopoint(dProp_smeared, dProp_stochEst, uProp_smeared, Base::proj_NO_PROJECTOR);
    C2_N_stochD += C2_N_stochD_tmp;
    C2_N_stochD *= 0.5;

    int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]};
    C2_N.prepareMomentumProjection(sourcePos);

    std::vector< Core::BaryonCorrelator > all_corrsP(C2_N.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_stochD(C2_N_stochD.momentumProjection(momenta));
    std::vector< Core::BaryonCorrelator > all_corrsP_stochU(C2_N_stochU.momentumProjection(momenta));

    std::ofstream *fout = NULL;
    if (weave.isRoot())
    {
      fout =  new std::ofstream("output_2point_neutron.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP[I].setOffset(timeslice_source);
        all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout->open("output_2point_neutron_stochD.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP_stochD[I].setOffset(timeslice_source);
        all_corrsP_stochD[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      fout->open("output_2point_neutron_stochU.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsP_stochU[I].setOffset(timeslice_source);
        all_corrsP_stochU[I].printWithMomentum_full(*fout, momenta[I]);
      }
      fout->close();
      delete fout;
    }

    C2_N_stochD.setOffset(timeslice_source);
    C2_N_stochD *= Base::proj_PARITY_PLUS_STD;
    if (weave.isRoot())
    {
      fout = new std::ofstream("output_2point_neutron_stochD_projected.dat");
      *fout << C2_N_stochD << std::endl;
      fout->close();
      delete fout;
    }
    C2_N_stochU.setOffset(timeslice_source);
    C2_N_stochU *= Base::proj_PARITY_PLUS_STD;
    if (weave.isRoot())
    {
      fout = new std::ofstream("output_2point_neutron_stochU_projected.dat");
      *fout << C2_N_stochU << std::endl;
      fout->close();
      delete fout;
    }
    C2_N.setOffset(timeslice_source);
    C2_N *= Base::proj_PARITY_PLUS_STD;
    if (weave.isRoot())
    {
      fout = new std::ofstream("output_2point_neutron_projected.dat");
      *fout << C2_N << std::endl;
      fout->close();
      delete fout;
    }
    weave.barrier();
    time (&end);
    dif = difftime (end,start);
    if (weave.isRoot())
      std::cout << "done in " << dif << " seconds\n" << std::endl;
  }

  weave.barrier();
  if (weave.isRoot())
    std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  return EXIT_SUCCESS;
}
