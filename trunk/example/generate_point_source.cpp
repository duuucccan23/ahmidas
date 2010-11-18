/*
   This piece of code generates a (smeared) point source
*/

//Ahmidas
#include <L0/Ahmidas.h>

// C++ complex library
#include <complex>

// C++ string library
#include <cstring>

// this we need for storing the input we read from an input file
#include <map>
#include <vector>

// this we need for C++ style output
#include <iomanip>
#include <iostream>


/* *** ahmidas interfaces *** */

// representation of Dirac gamma matrices
#include <L0/Dirac/Gamma.h>

// lattice structure
#include <L0/Core/Field.h>

// propagator structure
#include <L0/Core/Propagator.h>

// IO interface
#include <L1/Tool/IO.h>

// input file reader interface
#include <L2/Input/FileReader.h>


// interfaces needed for smearing
#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>


int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);


  // lattice size, to be read from input file (and therefore initialized to zero)
  size_t L_tmp(0);
  size_t T_tmp(0);


  /* those are the containers to be filled by the input file reader */

  // this one contains floating point variables like kappa, mu, ...
  std::map< std::string, double > floats;
  // this one contains the file names (themselves stored in containers) of propagators (or sources, gauge fields, ...)
  std::vector< std::vector< std::string > > files;

  std::vector< size_t * > positions;
  std::vector< int > operators;

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


  if(weave.isRoot())
  {
    std::cout << "\nThe following files are going to be created:" << std::endl;

    // there should only be one container in files, which can be accessed by files[0]
    // (similar to accessing an object in a C array)
    for (size_t fileIndex=0; fileIndex<files[0].size(); fileIndex++)
    {
      std::cout << (files[0])[fileIndex] << std::endl;
    }
  }


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


  Core::Field< QCD::Gauge > *gauge_field;

  if (Jac_iterations > 0)
  {
    gauge_field = new Core::Field< QCD::Gauge > (L, T);
    if (weave.isRoot())
      std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
    Tool::IO::load(gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
    if (weave.isRoot())
      std::cout << "done.\n" << std::endl;

    // smear the gauge field
    Smear::APE APE_tool(APE_alpha);
    APE_tool.smear(*gauge_field, APE_iterations);
    if (weave.isRoot())
      std::cout << "gauge field smeared successfully\n" << std::endl;
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
    std::cout << std::scientific << "norm of point source: " << norm << std::endl;
    if(norm != norm) // happens if norm is 'nan'
      exit(1);
  }

  weave.barrier();

  Tool::IO::save(&point_source, pointSourceFiles, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "point source saved successfully\n" << std::endl;


  // leave main function
  return EXIT_SUCCESS;
}
