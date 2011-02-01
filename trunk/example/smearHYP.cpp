
#include <string>
#include <iomanip>
#include <iostream>

#include <L0/Ahmidas.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Smear/HYP.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>

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


  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */

  // create input file reader, the name of the input file has to be passed as a parameter
  Input::FileReader reader("./smearHYP_input.xml");

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


  double const HYP_alpha1     = floats["HYP_param1"];
  double const HYP_alpha2     = floats["HYP_param2"];
  double const HYP_alpha3     = floats["HYP_param3"];
  size_t const HYP_iterations = size_t(floats["HYP_steps"]);

  if (weave.isRoot())
  {
    std::cout << "HYP smearing: parameters = " << HYP_alpha1 << ", " << HYP_alpha2 << ", " << HYP_alpha3
              << ", iterations = " << HYP_iterations << std::endl;
  }

  // for simpler use let us reference the file names via convenient names
  std::string const gaugeFieldInputFile((files[0])[0]);
  std::string const gaugeFieldOutputFile((files[1])[0]);

  if (weave.isRoot())
  {
    std::cout << "will read gauge field configuration from " << gaugeFieldInputFile
              << "\nand write smeared configuration to    "  << gaugeFieldOutputFile <<std::endl;
  }

  Core::Field< QCD::Gauge > gauge_field(L, T);
  Tool::IO::load(&gauge_field, gaugeFieldInputFile, Tool::IO::fileILDG);

  std::cout << "Showing the gauge link in the x-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_X];
  std::cout << "Showing the gauge link in the y-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_Y];
  std::cout << "Showing the gauge link in the z-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_Z];
  std::cout << "Showing the gauge link in the t-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_T];


  Smear::HYP hyp = Smear::HYP(HYP_alpha1, HYP_alpha2, HYP_alpha3);

  hyp.smear(gauge_field, HYP_iterations);

  std::cout << "\nShowing the smeared gauge link in the x-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_X];
  std::cout << "\nShowing the smeared gauge link in the y-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_Y];
  std::cout << "\nShowing the smeared gauge link in the z-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_Z];
  std::cout << "\nShowing the smeared gauge link in the t-up direction from point (0,0,0,0).\n";
  std::cout << gauge_field[0][Base::idx_T];

  Tool::IO::save(&gauge_field, gaugeFieldOutputFile, Tool::IO::fileILDG);

  return EXIT_SUCCESS;
}
