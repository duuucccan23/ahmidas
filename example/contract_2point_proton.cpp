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
#define __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__


int main(int argc, char **argv)
{
  Ahmidas ahmidas(&argc, &argv);

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./contract_2point_proton_input.xml");

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

//   Core::BaryonCorrelator C2_P_tm = Contract::proton_twopoint(forwardProp_u, forwardProp_u, forwardProp_d, Base::proj_PARITY_PLUS_TM);
//   C2_P_tm.deleteField();
//   C2_P_tm.setOffset(timeslice_source);
//   if (weave.isRoot())
//   {
//     std::ofstream fout("output_2point.dat");
//     fout << C2_P_tm;
//     fout.close();
//   }

  forwardProp_u.rotateToPhysicalBasis(true);
  forwardProp_d.rotateToPhysicalBasis(false);

  Core::BaryonCorrelator C2_P = Contract::proton_twopoint(forwardProp_u, forwardProp_u, forwardProp_d, Base::proj_NO_PROJECTOR);

  C2_P.deleteField();
  // C2_P.rotateToPhysicalBasis();

  C2_P.setOffset(timeslice_source);
  if (weave.isRoot())
  {
    std::ofstream fout("output_2point_unprojected.dat");
    for(size_t t = 0; t < T; t++)
    {
      fout.precision(8);
      fout.width(3);
      fout << t << "  " << std::scientific << std::showpos
           << C2_P[t][ 0].real() << "  " << C2_P[t][ 0].imag() << "  "
           << C2_P[t][ 1].real() << "  " << C2_P[t][ 1].imag() << "  "
           << C2_P[t][ 2].real() << "  " << C2_P[t][ 2].imag() << "  "
           << C2_P[t][ 3].real() << "  " << C2_P[t][ 3].imag() << "  "
           << std::endl;
      fout.width(3);
      fout << t << "  " << std::scientific << std::showpos
           << C2_P[t][ 4].real() << "  " << C2_P[t][ 4].imag() << "  "
           << C2_P[t][ 5].real() << "  " << C2_P[t][ 5].imag() << "  "
           << C2_P[t][ 6].real() << "  " << C2_P[t][ 6].imag() << "  "
           << C2_P[t][ 7].real() << "  " << C2_P[t][ 7].imag() << "  "
           << std::endl;
      fout.width(3);
      fout << t << "  " << std::scientific << std::showpos
           << C2_P[t][ 8].real() << "  " << C2_P[t][ 8].imag() << "  "
           << C2_P[t][ 9].real() << "  " << C2_P[t][ 9].imag() << "  "
           << C2_P[t][10].real() << "  " << C2_P[t][10].imag() << "  "
           << C2_P[t][11].real() << "  " << C2_P[t][11].imag() << "  "
           << std::endl;
      fout.width(3);
      fout << t << "  " << std::scientific << std::showpos
           << C2_P[t][12].real() << "  " << C2_P[t][12].imag() << "  "
           << C2_P[t][13].real() << "  " << C2_P[t][13].imag() << "  "
           << C2_P[t][14].real() << "  " << C2_P[t][14].imag() << "  "
           << C2_P[t][15].real() << "  " << C2_P[t][15].imag() << "  "
           << std::endl;

    }
    fout.close();
  }

  if (weave.isRoot())
    std::cout << "contractions performed and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
