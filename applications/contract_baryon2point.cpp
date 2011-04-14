// $Id: contract_nucleon3point.cpp 616 2011-04-07 10:56:18Z dinter@ifh.de $



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

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);
  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./contract_baryon2point_input.xml");

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



  size_t const * const source_position = positions[0];
  size_t const timeslice_source = source_position[Base::idx_T] % T;
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;

  // make sure the boundary is not crossed by source-sink correlaton function
   size_t const timeslice_boundary = (timeslice_source + (T/2)) % T;
//  size_t const timeslice_boundary = (timeslice_sink + 2) % T;
  if (weave.isRoot())
    std::cout << "timeslice (boundary) = " << timeslice_boundary << std::endl;

  std::vector< std::string > const &propfilesD(files[0]);
  std::vector< std::string > const &propfilesU(files[1]);
  std::vector< std::string > const &gaugeFieldFiles(files[2]);

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
         all_corrsP[I].printWithMomentum_full(*fout, momenta[I]);
        //all_corrsP[I] *= Base::proj_PARITY_PLUS_STD;
        //all_corrsP[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
      fout = new std::ofstream("output_2point_neutron.dat");
      for(size_t I=0; I<momenta.size(); I++)
      {
        all_corrsN[I].setOffset(timeslice_source);
        all_corrsN[I].printWithMomentum_full(*fout, momenta[I]);
        //all_corrsN[I] *= Base::proj_PARITY_PLUS_STD;
        //all_corrsN[I].printWithMomentum(*fout, momenta[I]);
      }
      fout->close();
    }

    for(size_t I=0; I<momenta.size(); I++)
      delete [] momenta[I];
    momenta.clear();
    delete fout;

    weave.barrier();

  }


  if (weave.isRoot())
    std::cout << "contractions performed and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
