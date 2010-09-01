
#include <complex>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>

#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Propagator.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>
#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>


int main(int argc, char **argv)
{


  // lattice size, to be read from input file (and therefore initialized to zero)
  size_t L_tmp(0);
  size_t T_tmp(0);


  /* those are the containers to be filled by the input file reader */

  // this one contains floating point variables like kappa, mu, ...
  std::map< std::string, double > floats;
  // this one contains the file names (themselves stored in containers) of propagators (or sources, gauge fields, ...)
  std::vector< std::vector< std::string > > files;

  std::vector< size_t * > positions;
  std::map< std::string, int > operators;

  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */

  // create input file reader, the name of the input file has to be passed as a parameter
  Input::FileReader reader("./generate_point_source_input.xml");

  // get input parameters
  // note: this is how to invoke a member function of an object in C++:
  // <name of object>.<name of function>(<parameter list>)
  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);

  const size_t L(L_tmp);
  const size_t T(T_tmp);

  /* ****************************************** */

  // this is needed if we want to have the output (i.e. to the standard output) done by only
  // one process in the parallel version
  Base::Weave weave(L, T);

  // that's how writing to the standard output works in C++
  // note: the "<<" operator also works for most of the ahmidas objects like SU3::Spinor or QCD::Tensor
  if(weave.isRoot())
    std::cout << "\nLattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  // that's how one can access the values in the map "floats"
  double kappa = floats["kappa"];
  double mu    = floats["mu"];

  if(weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

  size_t const * const source_position = positions[0];

  if(weave.isRoot())
    std::cout << "\nsource position: " << source_position[0] << " " << source_position[1] << " "
                                       << source_position[2] << " " << source_position[3] << std::endl;



  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = size_t(floats["APE_steps"]);
  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = size_t(floats["Jac_steps"]);

  if (weave.isRoot())
  {
    std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
    std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
  }

  // for simpler use let us reference the file names via convenient names
  std::vector< std::string > const & pointSourceFiles(files[0]);
  std::vector< std::string > const & gaugeFieldFiles(files[1]);

  double const TIME_FACTOR(1.0/CLOCKS_PER_SEC);
  double time_read;
  double time_smear;
  double time_write;


  Core::Field< QCD::Gauge > *gauge_field;



  if (Jac_iterations > 0)
  {

    gauge_field = new Core::Field< QCD::Gauge > (L, T);
    if (weave.isRoot())
      std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... " << std::flush;

    weave.barrier();
    double const time_1a = double(clock());

    Tool::IO::load(gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);

    weave.barrier();
    double const time_1b = double(clock());
    time_read = (time_1b - time_1a) * TIME_FACTOR;

    if (weave.isRoot())
    std::cout << "done in "
              << std::scientific << time_read
              << " seconds.\n" << std::endl;


    Smear::APE APE_tool(APE_alpha);

    weave.barrier();
    double const time_2a(clock());

    // smear the gauge field
    APE_tool.smear(*gauge_field, APE_iterations);

    weave.barrier();
    double const time_2b = double(clock());

    time_smear = (time_2b - time_2a) * TIME_FACTOR;

    if (weave.isRoot())
      std::cout << "gauge field smeared successfully in "
                << std::scientific << time_smear
                << " seconds\n" << std::endl;
  }

  Core::Propagator point_source(L, T);
  point_source *= 0.0;

  size_t localIndex = weave.globalCoordToLocalIndex(source_position[Base::idx_X],
                                                    source_position[Base::idx_Y],
                                                    source_position[Base::idx_Z],
                                                    source_position[Base::idx_T]);
  if (localIndex != weave.localVolume())
  {
    point_source[localIndex] = QCD::Tensor::identity();
  }
  if (Jac_iterations > 0)
  {
    point_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field);
    delete gauge_field;
    if(weave.isRoot())
      std::cout << "point source smeared successfully\n" << std::endl;
  }
  // to prevent NANs in optimized code on timeslices where everything should be zero
  point_source.select_timeslice(source_position[Base::idx_T]);
  double norm(point_source.norm());
  if (weave.isRoot())
  {
    std::cout.precision(10);
    std::cout << std::scientific << "norm of point source: " << norm << std::endl << std::endl;
    if(norm != norm) // happens if norm is 'nan'
      exit(1);
  }


  weave.barrier();
  double const time_3a = double(clock());


  Tool::IO::save(&point_source, pointSourceFiles, Tool::IO::fileSCIDAC);

  weave.barrier();
  double const time_3b = double(clock());

  time_write = (time_3b - time_3a) * TIME_FACTOR;

  if (weave.isRoot())
    std::cout << "point source saved successfully in "
              << std::scientific << time_write
              << " seconds\n" << std::endl;

  if (weave.isRoot())
  {
    std::cerr.precision(2);
    std::cerr << std::scientific << time_read  << "  "
              << std::scientific << time_smear << "  "
              << std::scientific << time_write << std::endl;
  }

  // leave main function
  return EXIT_SUCCESS;
}
