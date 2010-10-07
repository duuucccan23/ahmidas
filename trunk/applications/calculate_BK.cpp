/*
  This code calculates all three-point and two-point functions to measure B_K
  It uses the one-end-trick, assuming that the propagators are stochastic timeslice propagators
*/

// NOTE:
// the multi-mass solver of the tmLQCD package writes D^{dagger}D to the files.
// the isospin components have to be reconstructed 
// by application of D or D^{dagger}, respectively.

// personal note: replace Core::Propagator by Core::StochasticPropagator<4> in the future

#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Print.h>
#include <L0/Ahmidas.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>


int main(int argc, char **argv)
{

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./calculate_BK_input.xml");

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
  double mu_d    = floats["mu_d"];
  double mu_s    = floats["mu_s"];
  if (weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu_d = " << mu_d << ", mu_s = " << mu_s << std::endl;

  size_t const t_L = size_t(floats["t_L"]);
  assert(t_L >= 0 && t_L < T);
  size_t const t_R = size_t(floats["t_R"]);
  assert(t_R >= 0 && t_R < T);

  if (weave.isRoot())
    std::cout << "timeslice (left)  = " << t_L << std::endl;
  if (weave.isRoot())
    std::cout << "timeslice (right) = " << t_R << std::endl;

  std::vector< std::string > const &propfilesD_L(files[0]);
  std::vector< std::string > const &propfilesS_L(files[1]);
  std::vector< std::string > const &propfilesD_R(files[2]);
  std::vector< std::string > const &propfilesS_R(files[3]);
  std::vector< std::string > const &gaugeFieldFiles(files[4]);



/* ******************************************************************************* */
/* ******** load gauge field ***************************************************** */
/* ******************************************************************************* */


  Core::Field< QCD::Gauge > gauge_field(L, T);
  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;


/* ******************************************************************************* */
/* ******** load and prepare propagators from timeslice t_L ********************** */
/* ******************************************************************************* */

  // temporary Propagator just just for loading the MMSolver output
  Core::Propagator *tmpProp = new Core::Propagator(L, T);

  Tool::IO::load(tmpProp, propfilesD_L, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "1st propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp_plus_L(L, T);
  Core::Propagator dProp_minus_L(L, T);
  // reconstruct the full doublet from the combined MMS output, being ~ (DD^dagger)^-1
  tmpProp->reconstruct_doublet(dProp_plus_L, dProp_minus_L, gauge_field, kappa, mu_d);

  Tool::IO::load(tmpProp, propfilesS_L, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "2nd propagator successfully loaded\n" << std::endl;

  Core::Propagator sProp_plus_L(L, T);
  Core::Propagator sProp_minus_L(L, T);
  tmpProp->reconstruct_doublet(sProp_plus_L, sProp_minus_L, gauge_field, kappa, mu_s);

/* ******************************************************************************* */
/* ******** load and prepare propagators from timeslice t_R ********************** */
/* ******************************************************************************* */


  Tool::IO::load(tmpProp, propfilesD_R, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
      std::cout << "3rd propagator successfully loaded\n" << std::endl;

  Core::Propagator dProp_plus_R(L, T);
  Core::Propagator dProp_minus_R(L, T);
  tmpProp->reconstruct_doublet(dProp_plus_R, dProp_minus_R, gauge_field, kappa, mu_d);

  Tool::IO::load(tmpProp, propfilesS_R, Tool::IO::fileSCIDAC);
  if (weave.isRoot())
    std::cout << "4th propagator successfully loaded\n" << std::endl;

  Core::Propagator sProp_plus_R(L, T);
  Core::Propagator sProp_minus_R(L, T);
  tmpProp->reconstruct_doublet(sProp_plus_R, sProp_minus_R, gauge_field, kappa, mu_s);

  delete tmpProp;


// notes: construct objects O as in notes: multiply Dirac structure, dagger, and Tensor multiplication involved

  //gamma structure is: Gamma = gamma_mu (1 - gamma_5)

  Dirac::Gamma< 4 > gamma4;
  Dirac::Gamma< 54 > gamma5gamma4;

  // now we custruct some intermediate objects (called O_i(x,t) in the notes)

  // note: worry about efficiency when it works, maybe first taking 
  // the dagger of the d propagator might be advantageous, 
  // since we actually never need the undaggered propagator
  Core::Propagator O_R = gamma4 * dProp_plus_R;
  O_R += gamma5gamma4 * dProp_plus_R;
  O_R.dagger();
  O_R.rightMultiply(sProp_plus_R); // actually, is it right or left here? check again!


  Core::Propagator O_L = gamma4 * dProp_minus_L;
  O_L += gamma5gamma4 * dProp_minus_L;
  O_L.dagger();
  O_L.rightMultiply(sProp_plus_L); // actually, is it right or left here? check again!

  // there are more combinations (always with one 'minus' and three 'plus' propagators),
  // which can be written down straightforwardly

  // the disconnected part is just the product of the individual traces
  Core::Field < std::complex < double > > correlator = O_R.trace();
  correlator *= O_L.trace();
  correlator *= -1.0; // note that the sign is negative

  // now we can (matrix-)multiply O_L and O_R and take the trace thereof, being the connected part
  O_R.leftMultiply(O_L);
  correlator += O_R.trace();
  
  // this is not implemented yet ...

//  Core::Correlator< std::complex< double > > my_correlator(correlator); 
//  my_correlator.sumOverTimeslices();

//  std::ofstream fout(BK_output_3point.dat);
  if (weave.isRoot())
  {
//    fout << my_correlator;
//    fout.close();
  }

//  we also need a twopoint function (do it when constructing O_L and O_R, but just use gamma5 instead of gamma structure)

  if (weave.isRoot())
    std::cout << "contractions performed successfully\n" << std::endl;
  
  return EXIT_SUCCESS;
}
