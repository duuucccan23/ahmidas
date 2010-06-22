/*
   This piece of code generates stochastic timeslice sources with Z(4) noise
*/

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

// needed for initialization of random number generator
#include <L0/Base/Random.h>
#include <L0/Base/Z2.h>

// interfaces needed for smearing
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


  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */

  // create input file reader, the name of the input file has to be passed as a parameter
  Input::FileReader reader("../example/generate_source_input.xml");

  // get input parameters
  // note: this is how to invoke a member function of an object in C++:
  // <name of object>.<name of function>(<parameter list>)
  reader.initializeParameters(L_tmp, T_tmp, files, floats);

  const size_t L(L_tmp);
  const size_t T(T_tmp);

  /* ****************************************** */

  // this is needed if we want to have the output (i.e. to the standard output) done by only
  // one process in the parallel version
  Base::Weave weave(L, T);

  // ############################################################################################
  // ######### random number generation depend on the number of cores used ######################
  // ######### safe though for timeslice sources if parallelization only in (eucl.) time ########
  // ######### you can comment out the followind lines if you don't care ########################
  // ############################################################################################

  if(weave.isRoot())
    std::cerr << "FIX THIS: random number generation is safe only for scalar code so far,\n";
    std::cerr <<  "result would depend on the number of cores used!" << std::endl;
  weave.barrier();
  assert(weave.isRoot()); // this will fail if there's more than one process!

  // ##########################################################################################
  // ##########################################################################################


  // that's how writing to the standard output works in C++
  // note: the "<<" operator also works for most of the ahmidas objects like SU3::Spinor or QCD::Tensor
  if(weave.isRoot())
    std::cout << "\nLattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  // that's how one can access the values in the map "floats"
  double kappa = floats["kappa"];
  double mu    = floats["mu"];

  if(weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

  size_t t_src = size_t(floats["timesliceSource"]); 
  if(weave.isRoot())
    std::cout << "\nsource timeslice: " << t_src << std::endl;



  uint64_t const rSeed = uint64_t(floats["seed"]);
  if(weave.isRoot())
    std::cout << "\nrandom seed: " << rSeed << std::endl;
  Base::Z2::instance(1.0, rSeed);

  if (rSeed > 0)
  {
    double *tmp = new double[512];
    std::generate_n(tmp, 512, Base::Random::Z2);
    delete [] tmp;
  }

  if(weave.isRoot())
    std::cout << "\nThe following files are going to be created:" << std::endl;

  // there should only be one container in files, which can be accessed by files[0]
  // (similar to accessing an object in a C array)
  for (size_t fileIndex=0; fileIndex<files[0].size(); fileIndex++)
  {
    std::cout << (files[0])[fileIndex] << std::endl;
  }

  Base::SourcePolarization pol_tmp;
  if(int(floats["SourcePolarization"]) == Base::sou_FULLY_POLARIZED)
  {
    pol_tmp = Base::sou_FULLY_POLARIZED;
  }
  else if(int(floats["SourcePolarization"]) == Base::sou_PARTLY_POLARIZED)
  {
    pol_tmp = Base::sou_PARTLY_POLARIZED;
  }
  else
  {
    if(weave.isRoot())
      std::cerr << "source polarization " << floats["SourcePolarization"] << " unknown or not implemented" << std::endl;
    exit(1);
  }
  Base::SourcePolarization const polarization(pol_tmp);

  Base::SourceColorState col_tmp;
  if(int(floats["SourceColorState"]) == Base::sou_PURE)
  {
   col_tmp = Base::sou_PURE;
  }
  else if(int(floats["SourceColorState"]) == Base::sou_GENERIC)
  {
   col_tmp = Base::sou_GENERIC;
  }
  else
  {
    if(weave.isRoot())
      std::cerr << "source color state " << floats["SourceColorState"] << " unknown or not implemented" << std::endl;
    exit(1);
  }
  Base::SourceColorState const colorState(col_tmp);

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
  std::vector< std::string > const & stochasticSourceFiles(files[0]);
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
    APE_tool.smear(*gauge_field, APE_iterations, t_src);
    if (weave.isRoot())
      std::cout << "gauge field smeared successfully\n" << std::endl;
  }

  // version 1: spin (Dirac) and color dilution
  if(polarization == Base::sou_FULLY_POLARIZED && colorState == Base::sou_PURE)
  {
    Core::StochasticSource< 12 > stochastic_source(L, T, polarization, colorState, t_src);

    if (Jac_iterations > 0)
    {
      stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field, t_src);
      delete gauge_field;
      if(weave.isRoot())
        std::cout << "stochastic source smeared successfully\n" << std::endl;
    }
    Tool::IO::save(&stochastic_source, stochasticSourceFiles, Tool::IO::fileSCIDAC);
    if (weave.isRoot())
      std::cout << "stochastic source saved successfully\n" << std::endl;
  }
  // version 2: spin (Dirac) dilution only
  else if(polarization == Base::sou_FULLY_POLARIZED && colorState == Base::sou_GENERIC)
  {
    Core::StochasticSource< 4 > stochastic_source(L, T, polarization, colorState, t_src);

    if (Jac_iterations > 0)
    {
      stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field, t_src);
      delete gauge_field;
      if(weave.isRoot())
        std::cout << "stochastic source smeared successfully\n" << std::endl;
    }
    Tool::IO::save(reinterpret_cast< Core::StochasticPropagator< 4 > * >(&stochastic_source), stochasticSourceFiles, Tool::IO::fileSCIDAC);
    if (weave.isRoot())
      std::cout << "stochastic source saved successfully\n" << std::endl;
  }
  // version 2: spin (Dirac) dilution only
  else if(polarization == Base::sou_PARTLY_POLARIZED && colorState == Base::sou_GENERIC)
  {
    Core::StochasticSource< 1 > stochastic_source(L, T, polarization, colorState, t_src);

    if (Jac_iterations > 0)
    {
      stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field, t_src);
      delete gauge_field;
      if(weave.isRoot())
        std::cout << "stochastic source smeared successfully\n" << std::endl;
    }
    Tool::IO::save(reinterpret_cast< Core::StochasticPropagator< 1 > * >(&stochastic_source), stochasticSourceFiles, Tool::IO::fileSCIDAC);
    if (weave.isRoot())
      std::cout << "stochastic source saved successfully\n" << std::endl;
  }
  else
  {
    if(weave.isRoot())
      std::cerr << "source polarization and color state combination not implemented" << std::endl;
    exit(1);
  }


  // leave main function
  return EXIT_SUCCESS;
}
