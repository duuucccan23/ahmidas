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

// use operator involving gamma_1
//#define __GAMMA_1__
// use operator involving gamma_2
 //#define __GAMMA_2__
// use operator involving gamma_3
// #define __GAMMA_3__
// use operator involving gamma_4
// #define __GAMMA_4__

#define __O44__

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./contract_3point_sequential_alternative_input.xml");

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

  size_t const timeslice_source = size_t(floats["timesliceSource"]);
  if (weave.isRoot())
    std::cout << "timeslice (source) = " << timeslice_source << std::endl;

  size_t const t_op = size_t(floats["timesliceInsertion"]);
  if (weave.isRoot())
    std::cout << "timeslice (operator insertion) = " << t_op << std::endl;

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

  Core::Propagator uProp(L, T);
  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark forward propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp(L, T);
  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark forward propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_FW__
  dProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  uProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
#endif

  Core::Propagator seqProp_u(L, T);
  Core::Propagator seqProp_d(L, T);

  Tool::IO::load(&seqProp_u, seqPropFilesU, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "u quark backward propagator successfully loaded\n" << std::endl;

  Tool::IO::load(&seqProp_d, seqPropFilesD, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "d quark backward propagator successfully loaded\n" << std::endl;


#ifdef __COMPENSATE_UNIFORM_BOUNDARY_CONDITIONS_BW__
  seqProp_u.changeBoundaryConditions_uniformToFixed(t_op, timeslice_boundary);
  seqProp_d.changeBoundaryConditions_uniformToFixed(t_op, timeslice_boundary);
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

  uProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  dProp.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  seqProp_u.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
  seqProp_d.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

  if (weave.isRoot())
    std::cout << "propagators and gauge field smeared successfully\n" << std::endl;


  {

    if (weave.isRoot())
      std::cout << "\n calculating 3-point function(s) \n" << std::endl;

    Core::BaryonCorrelator C3p_d1 = Contract::proton_twopoint(uProp, uProp, seqProp_d, Base::proj_NO_PROJECTOR);
    Core::BaryonCorrelator C3p_u1 = Contract::proton_twopoint(uProp, seqProp_u, dProp, Base::proj_NO_PROJECTOR);
    Core::BaryonCorrelator C3p_u2 = Contract::proton_twopoint(seqProp_u, uProp, dProp, Base::proj_NO_PROJECTOR);

    C3p_u1.deleteField();
    C3p_d1.deleteField();
    C3p_u2.deleteField();

    C3p_u1 += C3p_u2;

#ifdef __GAMMA_1__
    Core::BaryonCorrelator C3p_u(C3p_u1);
    Core::BaryonCorrelator C3p_d(C3p_d1);
    C3p_u *= Base::proj_1_PLUS_TM;
    C3p_d *= Base::proj_1_PLUS_TM;
#endif
#ifdef __GAMMA_2__
    Core::BaryonCorrelator C3p_u(C3p_u1);
    Core::BaryonCorrelator C3p_d(C3p_d1);
    C3p_u *= Base::proj_2_PLUS_TM;
    C3p_d *= Base::proj_2_PLUS_TM;
#endif
#ifdef __GAMMA_3__
    Core::BaryonCorrelator C3p_u(C3p_u1);
    Core::BaryonCorrelator C3p_d(C3p_d1);
    C3p_u *= Base::proj_3_PLUS_TM;
    C3p_d *= Base::proj_3_PLUS_TM;
#endif
#ifdef __GAMMA_4__
    Core::BaryonCorrelator C3p_u(C3p_u1);
    Core::BaryonCorrelator C3p_d(C3p_d1);
    C3p_u *= Base::proj_PARITY_PLUS_TM;
    C3p_d *= Base::proj_PARITY_PLUS_TM;
#endif

#ifdef __O44__
    Core::BaryonCorrelator C3p_u(C3p_u1);
    Core::BaryonCorrelator C3p_d(C3p_d1);
    C3p_u *= Base::proj_PARITY_PLUS_TM;
    C3p_d *= Base::proj_PARITY_PLUS_TM;
#endif


    C3p_u.setOffset(timeslice_source);
    C3p_d.setOffset(timeslice_source);

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
     #ifdef __O44__
     std::ofstream fout("output_3point_O44_subtr_uu.dat");
     #endif

     fout << C3p_u << std::endl;
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
     #ifdef __O44__
     fout.open("output_3point_O44_subtr_dd.dat");
     #endif


     fout << C3p_d << std::endl;
     fout.close();
    }
  }


  {
    if (weave.isRoot())
      std::cout << "\n calculating 2-point function \n" << std::endl;

    Core::BaryonCorrelator C2_P = Contract::proton_twopoint(uProp, uProp, dProp, Base::proj_PARITY_PLUS_TM);
    C2_P.setOffset(timeslice_source);
    if (weave.isRoot())
    {
      std::ofstream fout("output_2point.dat");
      std::cout << "proton two point" << std::endl;
      std::cout << C2_P << std::endl;
      fout << C2_P << std::endl;
      fout.close();
    }
  }


  {
    if (weave.isRoot())
      std::cout << "\n calculating 2-point function (unprojected) \n" << std::endl;

    uProp.rotateToPhysicalBasis(true);
    dProp.rotateToPhysicalBasis(false);

    Core::BaryonCorrelator C2_P = Contract::proton_twopoint(uProp, uProp, dProp, Base::proj_NO_PROJECTOR);

    C2_P.deleteField();

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
  }
  
  if (weave.isRoot())
    std::cout << "contractions performed and saved successfully\n" << std::endl;

  return EXIT_SUCCESS;
}
